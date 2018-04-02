
#include <algorithm>
#include "LongReadProbe.h"
#include "Util.h"
#include "KmerFeature.h"
#include "KmerThreshold.h"

thread_local ProbeParameters LongReadProbe::m_params;

// Search seeds with [static/dynamic] kmers. Noted by KuanWeiLee 20171027
void LongReadProbe::searchSeedsWithHybridKmers(const std::string& readSeq, SeedFeature::SeedVector& seedVec)
{
	const size_t readSeqLen = readSeq.length();
	int staticSize = m_params.startKmerLen;
	if((int)readSeqLen < staticSize) return;
	
	int* attribute = new int[readSeqLen];
	getSeqAttribute(readSeq, attribute);
	if(m_params.Manual) std::fill_n(attribute, readSeqLen, m_params.mode);
	
	//Search seeds; slide through the read sequence with hybrid-kmers. Noted by KuanWeiLee
	//[init/curr]Pos indicate the initial/current position of the static-kmer.
	for(size_t initPos = 0; initPos < readSeqLen; initPos++)
	{
		int dynamicMode = attribute[initPos];
		staticSize += m_params.kmerOffset[dynamicMode];
		bool isSeed = false, isRepeat = false;
		KmerFeature dynamicKmer = KmerFeature::kmerRec[staticSize][initPos];
		int maxFixedMerFreq = dynamicKmer.getFreq();
		size_t seedStartPos = initPos;
		for(size_t currPos = initPos; currPos < readSeqLen; currPos++)
		{
			int staticMode = attribute[currPos];
			const KmerFeature& staticKmer = KmerFeature::kmerRec[staticSize][currPos];
			if(isSeed)
			{
				char b = readSeq[(currPos + staticSize - 1)];
				dynamicKmer.expand(b);
			}
			float dynamicThreshold = KmerThreshold::Instance().get(dynamicMode, dynamicKmer.getSize());
			float staticThreshold  = KmerThreshold::Instance().get(staticMode,  staticKmer.getSize());
			float repeatThreshold  = staticThreshold * (5 - ((staticMode >> 1) << 2));
			bool isOnRepeat = (staticKmer.getFreq() >= repeatThreshold);
			float freqDiff = (float)staticKmer.getFreq()/maxFixedMerFreq;
			//Gerneral seed extension strategy.
			if	(
				   staticKmer.getFreq() < staticThreshold						//1.static frequency
				|| dynamicKmer.getFreq() < dynamicThreshold						//2.dynamic frequency(1)
				|| !dynamicKmer.isValid()										//2.dynamic frequency(2)
				|| dynamicKmer.getSize() > m_params.kmerLenUpBound				//3.over length
				)
			{
				if(isSeed && !staticKmer.getPseudo()) dynamicKmer.shrink(1);
				break;
			}
			//Kmer Hitchhike strategy.
			int isGiantRepeat = ((dynamicMode >> 1) & (staticMode >> 1)) + 1;
			if(isRepeat && freqDiff < (m_params.hhRatio/isGiantRepeat))			//4.hitchhiking kmer(1) (HIGH-->LOW)
			{
				dynamicKmer.shrink(1);
				initPos++;
				break;
			}
			else if(isOnRepeat && freqDiff > (isGiantRepeat/m_params.hhRatio))	//4.hitchhiking kmer(2) (LOW-->HIGH)
			{
				isSeed = false;
				initPos = currPos - 1;
				break;
			}
			isSeed = true;
			initPos = seedStartPos + dynamicKmer.getSize() - 1;
			isRepeat = isRepeat || isOnRepeat;
			maxFixedMerFreq = std::max(maxFixedMerFreq, staticKmer.getFreq());
		}
		//Low Complexity strategy.
		if(isSeed && !dynamicKmer.isLowComplexity())
		{
			SeedFeature newSeed(dynamicKmer.getWord(), seedStartPos, maxFixedMerFreq, isRepeat, staticSize, m_params.PBcoverage);
			newSeed.estimateBestKmerSize(m_params.indices);
			seedVec.push_back(newSeed);
		}
		staticSize -= m_params.kmerOffset[dynamicMode];
	}

	//Seed Hitchhike strategy.
	seedVec = removeHitchhikingSeeds(seedVec, attribute);
	if(m_params.DebugSeed)
	{
		std::ostream* pSeedWriter = createWriter(m_params.directory + "seed/" + m_params.readid + ".seed");
		write(*pSeedWriter, seedVec);
		delete pSeedWriter;
	}
	
	delete[] attribute;
}
//Sequence attribute is set dynamically using a sliding fixed-mer on each position of the sequence.
//Noted by KuanWeiLee 20180118
void LongReadProbe::getSeqAttribute(const std::string& seq, int* const attribute)
{
	const size_t seqLen = seq.length();
	std::fill_n(attribute, seqLen, 1);
	
	int range = 300;
	const int ksize = m_params.scanKmerLen;
	float repeatValue = KmerThreshold::Instance().get(2, ksize);
	
	int front = 0, fear = -1;
	int leftmost = (seqLen - 1), rightmost = 0, repeatcount = 0;
	std::map<int, int> set;
	
	for(size_t pos = 0; pos < seqLen; pos++)
	{
		int left = pos - (range >> 1);
		int right = pos + (range >> 1);
		left  = std::max(left, 0);
		right = std::min(right, (int)(seqLen - 1));
		while(fear < right)
		{
			fear++;
			KmerFeature* prev = nullptr;
			for(auto& iter : m_params.kmerPool)
			{
				KmerFeature::kmerRec[iter][fear] = KmerFeature(m_params.indices, seq, fear, iter, prev);
				prev = KmerFeature::kmerRec[iter].get() + fear;
			}
			const KmerFeature& inKmer = KmerFeature::kmerRec[ksize][fear];
			int freq = inKmer.isLowComplexity() ? 0 : inKmer.getFreq();
			int mode;
			if(freq == 0) mode = -1;
			else if(freq >= repeatValue) mode = 2;
			else mode = 1;
			set[mode]++;
		}
		while(front < left)
		{
			const KmerFeature& outKmer = KmerFeature::kmerRec[ksize][front];
			front++;
			int freq = outKmer.isLowComplexity() ? 0 : outKmer.getFreq();
			int mode;
			if(freq == 0) mode = -1;
			else if(freq >= repeatValue) mode = 2;
			else mode = 1;
			set[mode]--;
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
	
//	if((float)repeatcount/seqLen >= 0.5 && (float)(leftmost + (seqLen -rightmost))/seqLen <= 0.1)
//		std::fill_n(attribute, seqLen, 2);
}

//Kmer & Seed Hitchhike strategy would maitain seed-correctness, 
//once the sequence is stuck between the ambiguity from uniqu to repeatThreshold mode.
//Noted by KuanWeiLee 20180106
SeedFeature::SeedVector LongReadProbe::removeHitchhikingSeeds(SeedFeature::SeedVector initSeedVec, int const *attribute)
{
	if(initSeedVec.size() < 2) return initSeedVec;
	
	int ksize = m_params.startKmerLen + m_params.kmerOffset[2];
	float overFreq = KmerThreshold::Instance().get(2, ksize)*5;
	
	for(SeedFeature::SeedVector::iterator iterQuery = initSeedVec.begin(); (iterQuery + 1) != initSeedVec.end(); iterQuery++)
	{
		SeedFeature& query = *iterQuery;
		SeedFeature::SeedVector::iterator iterTarget = iterQuery + 1;
		//if(query.isHitchhiked) continue;
		int queryMode = attribute[query.seedStartPos];
		if(queryMode == 2 && query.maxFixedMerFreq >= overFreq) continue;
		
		for(; iterTarget != initSeedVec.end(); iterTarget++)
		{
			SeedFeature& target = *iterTarget;
			//if(target.isHitchhiked) continue;
			int	targetMode = attribute[target.seedStartPos];
		//	int isBoundary = (queryMode >> 1) ^ (targetMode >> 1);
		//	if((int)(target.seedStartPos - query.seedEndPos) > (m_params.repeatDis >> isBoundary)) break;
			if((int)(target.seedStartPos - query.seedEndPos) > m_params.repeatDis) break;
			if(targetMode == 2 && target.maxFixedMerFreq >= overFreq) continue;
			float freqDiff = (float)target.maxFixedMerFreq/query.maxFixedMerFreq;
			int isGiantRepeat = ((queryMode >> 1) & (targetMode >> 1)) + 1;
			
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