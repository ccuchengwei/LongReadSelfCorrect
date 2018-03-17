
#include <algorithm>
#include "LongReadProbe.h"
#include "Util.h"
#include "KmerFeature.h"
#include "KmerThreshold.h"

thread_local ProbeParameters LongReadProbe::m_params = ProbeParameters();

// Search seeds with [static/dynamic] kmers. Noted by KuanWeiLee 20171027
void LongReadProbe::searchSeedsWithHybridKmers(const std::string& readSeq, SeedFeature::SeedVector& seedVec)
{
	const size_t readSeqLen = readSeq.length();
	int staticKmerSize = m_params.startKmerLen;
	if((int)readSeqLen < staticKmerSize) return;
	
	int attribute[readSeqLen];
	getSeqAttribute(readSeq, attribute);
	if(m_params.Manual) std::fill_n(attribute, readSeqLen, m_params.mode);
	
	//Search seeds; slide through the read sequence with hybrid-kmers. Noted by KuanWeiLee
	//kmer[Start/Move]Pos indicate the starting/moving position of the current static-kmer.
	float** table = KmerThreshold::m_table;

	for(size_t kmerStartPos = 0; kmerStartPos < readSeqLen; kmerStartPos++)
	{
		int kmerStartType = attribute[kmerStartPos];
		staticKmerSize += m_params.kmerOffset[kmerStartType];
		bool isSeed = false, isRepeat = false;
		KmerFeature dynamicKmer = KmerFeature::kmerRec[staticKmerSize][kmerStartPos];
	//	KmerFeature dynamicKmer(&m_params.indices, readSeq, kmerStartPos, staticKmerSize);
		int maxFixedMerFreq = dynamicKmer.getFreq();
		size_t seedStartPos = kmerStartPos;
		for(size_t kmerMovePos = kmerStartPos; kmerMovePos < readSeqLen; kmerMovePos++)
		{
			int kmerMoveType = attribute[kmerMovePos];
			const KmerFeature& staticKmer = KmerFeature::kmerRec[staticKmerSize][kmerMovePos];
		//	const KmerFeature staticKmer(&m_params.indices, readSeq, kmerMovePos, staticKmerSize);
			if(isSeed)
			{
				char b = readSeq[(kmerMovePos + staticKmerSize - 1)];
				dynamicKmer.expand(b);
			}
			float dynamicThreshold = table[kmerStartType][dynamicKmer.getSize()];
		//	float dynamicThreshold = table[            1][dynamicKmer.getSize()]; //wait for delete
			float staticThreshold = table[kmerMoveType][staticKmerSize];
		//	float staticThreshold = table[           1][staticKmerSize]; //wait for delete
			float repeatThreshold = staticThreshold * (5 - ((kmerMoveType >> 1) << 2));
			bool isOnRepeat = (staticKmer.getFreq() >= repeatThreshold);
			float freqDiff = (float)staticKmer.getFreq()/maxFixedMerFreq;
			//Gerneral seed extension strategy.
			if	(
				   dynamicKmer.getSize() > m_params.kmerLenUpBound				//1.over length
				|| !dynamicKmer.isValid()										//2.dynamic frequency(1)
				|| dynamicKmer.getFreq() < dynamicThreshold						//2.dynamic frequency(2)
				|| staticKmer.getFreq() < staticThreshold						//3.static frequency
				)
			{
				if(isSeed && !staticKmer.getPseudo()) dynamicKmer.shrink(1);
				break;
			}
			//Kmer Hitchhike strategy.
			int isGiantRepeat = ((kmerStartType >> 1) & (kmerMoveType >> 1)) + 1;
			if(isRepeat && freqDiff < (m_params.hhRatio/isGiantRepeat))			//4.hitchhiking kmer(1) (HIGH-->LOW)
			{
				dynamicKmer.shrink(1);
				kmerStartPos++;
				break;
			}
			else if(isOnRepeat && freqDiff > (isGiantRepeat/m_params.hhRatio))	//4.hitchhiking kmer(2) (LOW-->HIGH)
			{
				isSeed = false;
				kmerStartPos = kmerMovePos - 1;
				break;
			}
			isSeed = true;
			kmerStartPos = seedStartPos + dynamicKmer.getSize() - 1;
			isRepeat = isRepeat || isOnRepeat;
			maxFixedMerFreq = std::max(maxFixedMerFreq, staticKmer.getFreq());
		}
		//Low Complexity strategy.
		if(isSeed && !dynamicKmer.isLowComplexity())
		{
			SeedFeature newSeed(dynamicKmer.getWord(), seedStartPos, maxFixedMerFreq, isRepeat, staticKmerSize, m_params.PBcoverage);
			newSeed.estimateBestKmerSize(m_params.indices);
			seedVec.push_back(newSeed);
		}
		staticKmerSize -= m_params.kmerOffset[kmerStartType];
	}

	//Seed Hitchhike strategy.
	seedVec = removeHitchhikingSeeds(seedVec, attribute);
	if(m_params.DebugSeed)
	{
		std::ostream* pSeedWriter = createWriter(m_params.directory + "seed/" + m_params.readid + ".seed");
		write(*pSeedWriter, seedVec);
		delete pSeedWriter;
	}
}
//Sequence attribute is set dynamically using a sliding fixed-mer on each position of the sequence.
//Noted by KuanWeiLee 20180118
void LongReadProbe::getSeqAttribute(const std::string& seq, int* const attribute)
{
	const size_t seqlen = seq.length();
	std::fill_n(attribute, seqlen, 1);
	
	int range = 300;
	int x = m_params.PBcoverage;
	int y = m_params.scanKmerLen;
	const int ksize = m_params.scanKmerLen;
//	float lowcov = KmerThreshold::calculate(0, x, y);
//	float unique = KmerThreshold::calculate(1, x, y);
	float repeat = KmerThreshold::calculate(2, x ,y);
	int front = 0, fear = -1;
	int leftmost = (seqlen - 1), rightmost = 0, repeatcount = 0;
	std::map<int, int> set;
	
	for(size_t pos = 0; pos < seqlen; pos++)
	{
		int left = pos - (range >> 1);
		int right = pos + (range >> 1);
		left = std::max(left, 0);
		right = std::min(right, (int)(seqlen - 1));
		while(fear < right)
		{
			fear++;
			KmerFeature* prev = nullptr;
			for(auto& iter : m_params.kmerPool)
			{
				KmerFeature::kmerRec[iter][fear] = KmerFeature(&m_params.indices, seq, fear, iter, prev);
				prev = KmerFeature::kmerRec[iter].get() + fear;
			}
			const KmerFeature& inKmer = KmerFeature::kmerRec[ksize][fear];
		//	Kmer inKmer(&m_params.indices, seq, ++fear, ksize);
			int freq = inKmer.isLowComplexity() ? 0 : inKmer.getFreq();
		//	int freq = inKmer.getFreq();
			int type;
			if(freq == 0) type = -1;
			else if(freq >= repeat) type = 2;
			else type = 1;
			set[type]++;
		}
		while(front < left)
		{
			const KmerFeature& outKmer = KmerFeature::kmerRec[ksize][front];
			front++;
		//	Kmer outKmer(&m_params.indices, seq, front++, ksize);
			int freq = outKmer.isLowComplexity() ? 0 : outKmer.getFreq();
		//	int freq = outKmer.getFreq();
			int type;
			if(freq == 0) type = -1;
			else if(freq >= repeat) type = 2;
			else type = 1;
			set[type]--;
		}
		int size = (right - left + 1) - set[-1];
		float ratio = (float)set[2]/size + 0.0005;
		if(ratio >= 0.02)
		{
			attribute[pos] = 2;
			repeatcount++;
			leftmost = std::min(leftmost, (int)pos);
			rightmost = std::max(rightmost, (int)pos);
		}
	}
	
	if((float)repeatcount/seqlen >= 0.5 && (float)(leftmost + (seqlen -rightmost))/seqlen <= 0.1)
		std::fill_n(attribute, seqlen, 2);
}

//Kmer & Seed Hitchhike strategy would maitain seed-correctness, 
//once the sequence is stuck between the ambiguity from uniqu to repeat mode.
//Noted by KuanWeiLee 20180106
SeedFeature::SeedVector LongReadProbe::removeHitchhikingSeeds(SeedFeature::SeedVector initSeedVec, int const *attribute)
{
	if(initSeedVec.size() < 2) return initSeedVec;
	
	int x = m_params.PBcoverage;
	int y = m_params.startKmerLen + m_params.kmerOffset[2];
	float overFreq = KmerThreshold::calculate(2, x ,y) * 5;
	
	for(SeedFeature::SeedVector::iterator iterQuery = initSeedVec.begin(); (iterQuery + 1) != initSeedVec.end(); iterQuery++)
	{
		SeedFeature& query = *iterQuery;
		SeedFeature::SeedVector::iterator iterTarget = iterQuery + 1;
		//if(query.isHitchhiked) continue;
		int queryType = attribute[query.seedStartPos];
		if(queryType == 2 && query.maxFixedMerFreq >= overFreq) continue;
		
		for(; iterTarget != initSeedVec.end(); iterTarget++)
		{
			SeedFeature& target = *iterTarget;
			//if(target.isHitchhiked) continue;
			int	targetType = attribute[target.seedStartPos];
		//	int isBoundary = (queryType >> 1) ^ (targetType >> 1);
		//	if((int)(target.seedStartPos - query.seedEndPos) > (m_params.repeatDis >> isBoundary)) break;
			if((int)(target.seedStartPos - query.seedEndPos) > m_params.repeatDis) break;
			if(targetType == 2 && target.maxFixedMerFreq >= overFreq) continue;
			float freqDiff = (float)target.maxFixedMerFreq/query.maxFixedMerFreq;
			int isGiantRepeat = ((queryType >> 1) & (targetType >> 1)) + 1;
			
			target.isHitchhiked = target.isHitchhiked || (query.isRepeat && freqDiff < (m_params.hhRatio/isGiantRepeat));	//HIGH --> LOW
			query.isHitchhiked = query.isHitchhiked || (target.isRepeat && freqDiff > (isGiantRepeat/m_params.hhRatio));	//LOW  --> HIGH
		}
	}
	
	SeedFeature::SeedVector finalSeedVec, outcastSeedVec;
	finalSeedVec.reserve(initSeedVec.size());
	outcastSeedVec.reserve(initSeedVec.size() >> 1);
	
	for(const auto& iter : initSeedVec)
	{
		if(iter.isHitchhiked)
			outcastSeedVec.push_back(iter);
		else
			finalSeedVec.push_back(iter);
	}
	
	if(m_params.DebugSeed)
	{
		std::ostream* pOutcastSeedWriter = createWriter(m_params.directory + "seed/error/" + m_params.readid + ".seed");
		write(*pOutcastSeedWriter, outcastSeedVec);
		delete pOutcastSeedWriter;
	}
	return finalSeedVec;
}

void LongReadProbe::write(std::ostream& outfile, const SeedFeature::SeedVector& seedVec)
{	
	for(const auto& iter : seedVec)
		outfile
		<< iter.seedStr << "\t"
		<< iter.maxFixedMerFreq << "\t" 
		<< iter.seedStartPos << "\t"
		<< (iter.isRepeat ? "Yes" : "No") << "\n";
}