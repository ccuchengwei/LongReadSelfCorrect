
#include <algorithm>
#include "LongReadProbe.h"
#include "Util.h"
#include "KmerFeature.h"
#include "KmerThreshold.h"


ProbeParameters::ProbeParameters(
		BWTIndexSet _indices,
		std::string _directory,
		int _startKmerLen,
		int _PBcoverage,
		int _mode,
		std::array<int, 3> _offset,
		std::set<int> _pool,
		bool _DebugSeed,
		bool _Manual)
:	indices(_indices),
	directory(_directory),
	startKmerLen(_startKmerLen),
	PBcoverage(_PBcoverage),
	mode(_mode),
	offset(_offset),
	pool(_pool),
	DebugSeed(_DebugSeed),
	Manual(_Manual){ }

ProbeParameters LongReadProbe::m_params;

thread_local std::string LongReadProbe::readid;

// Search seeds with [static/dynamic] kmers. Noted by KuanWeiLee 20171027
void LongReadProbe::searchSeedsWithHybridKmers(const std::string& readSeq, SeedFeature::SeedVector& seedVec)
{
	const size_t readSeqLen = readSeq.length();
	int staticSize = m_params.startKmerLen;
	if((int)readSeqLen < staticSize) return;
	
	int* attribute = new int[readSeqLen];
	getSeqAttribute(readSeq, attribute);
	if(m_params.Manual) std::fill_n(attribute, readSeqLen, m_params.mode);
	
	//Search seeds; slide through the sequence with hybrid-kmers. Noted by KuanWeiLee
	//[init/curr]Pos indicate the initial/current position of the static-kmer.
	for(size_t initPos = 0; initPos < readSeqLen; initPos++)
	{
		int dynamicMode = attribute[initPos];
		staticSize += m_params.offset[dynamicMode];
		KmerFeature dynamicKmer = KmerFeature::Log()[staticSize][initPos];
		bool isSeed = false, isRepeat = false;
		int maxFixedMerFreq = dynamicKmer.getFreq();
		size_t seedPos = initPos;
		for(size_t currPos = initPos; currPos < readSeqLen; currPos++)
		{
			int staticMode = attribute[currPos];
			const KmerFeature& staticKmer = KmerFeature::Log()[staticSize][currPos];
			if(staticKmer.isFake()) break;
			if(isSeed)
			{
				char b = readSeq[(currPos + staticSize - 1)];
				dynamicKmer.expand(b);
			}
			float dynamicThreshold = KmerThreshold::Instance().get(dynamicMode, dynamicKmer.getSize());
			float staticThreshold  = KmerThreshold::Instance().get(staticMode, staticKmer.getSize());
			float repeatThreshold  = (5 - ((staticMode >> 1) << 2))*staticThreshold;
			//Gerneral seed extension strategy.
			if	(
				   staticKmer.getFreq() < staticThreshold						//1.static frequency
				|| dynamicKmer.getFreq() < dynamicThreshold						//2.dynamic frequency(1)
				|| !dynamicKmer.isValid()										//2.dynamic frequency(2)
				|| dynamicKmer.getSize() > m_params.kmerLenUpBound				//3.over length
				)
			{
				if(isSeed) dynamicKmer.shrink(1);
				break;
			}
			//Kmer Hitchhike strategy.
			float freqDiff = (float)staticKmer.getFreq()/maxFixedMerFreq;
			if(freqDiff < m_params.hhRatio)						//4.hitchhiking kmer(1) (HIGH-->LOW)
			{
				initPos++;
				dynamicKmer.shrink(1);
				break;
			}
			else if(freqDiff > 1/m_params.hhRatio)				//4.hitchhiking kmer(2) (LOW-->HIGH)
			{
				initPos = currPos - 1;
				isSeed = false;
				break;
			}
			initPos = seedPos + dynamicKmer.getSize() - 1;
			isSeed = true;
			isRepeat |= (staticKmer.getFreq() >= repeatThreshold);
			maxFixedMerFreq = std::max(maxFixedMerFreq, staticKmer.getFreq());
		}
		//Low Complexity strategy.
		if(isSeed && !dynamicKmer.isLowComplexity())
		{
			seedVec.push_back(SeedFeature(dynamicKmer.getWord(), seedPos, maxFixedMerFreq, isRepeat, staticSize, m_params.PBcoverage));
			seedVec.back().estimateBestKmerSize(m_params.indices);
		}
		staticSize -= m_params.offset[dynamicMode];
	}

	//Seed Hitchhike strategy.
	seedVec = removeHitchhikingSeeds(seedVec);
	
	if(m_params.DebugSeed)
	{
		std::ostream* pSeedWriter = createWriter(m_params.directory + "seed/" + readid + ".seed");
		*pSeedWriter << seedVec;
		delete pSeedWriter;
	}
	
	delete[] attribute;
}
//Sequence attribute is set dynamically using a sliding fixed-mer on each position of the sequence.
//Noted by KuanWeiLee 20180118
void LongReadProbe::getSeqAttribute(const std::string& seq, int* const attribute)
{
	std::ostream* pAutoWriter = nullptr;
	if(m_params.DebugSeed)
		pAutoWriter = createWriter(m_params.directory + "extend/" + readid + ".log");
	const size_t seqLen = seq.length();
	std::fill_n(attribute, seqLen, 1);
	
	int range = 300;
	const int ksize = m_params.scanKmerLen;
	float repeatValue = KmerThreshold::Instance().get(2, ksize);
	
	int front = 0, fear = -1;
	int leftmost = (seqLen - 1), rightmost = 0;
	std::map<int, int> box; //-1 -> garbage; 0 -> lowcov(disable); 1 -> unique; 2 -> repeat
	
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
			for(auto& iter : m_params.pool)
			{
				KmerFeature::Log()[iter][fear] = KmerFeature(m_params.indices, seq, fear, iter, prev);
				prev = KmerFeature::Log()[iter].get() + fear;
			}
			const KmerFeature& inKmer = KmerFeature::Log()[ksize][fear];
			int freq = inKmer.isLowComplexity() ? -1 : inKmer.getFreq();
			int mode;
			if(freq < 0) mode = -1;
			else if(freq >= repeatValue) mode = 2;
			else mode = 1;
			box[mode]++;
		}
		while(front < left)
		{
			const KmerFeature& outKmer = KmerFeature::Log()[ksize][front];
			front++;
			int freq = outKmer.isLowComplexity() ? -1 : outKmer.getFreq();
			int mode;
			if(freq <= 0) mode = -1;
			else if(freq >= repeatValue) mode = 2;
			else mode = 1;
			box[mode]--;
		}
		int size = (right - left + 1) - box[-1];
		float ratio = (float)box[2]/size + 0.0005;
		if(m_params.DebugSeed)
			*pAutoWriter << pos << '\t' << ratio << '\n';
		if(ratio >= 0.02)
		{
			attribute[pos] = 2;
			leftmost = std::min(leftmost, (int)pos);
			rightmost = std::max(rightmost, (int)pos);
		}
	}
	delete pAutoWriter;
}

//Kmer & Seed Hitchhike strategy would maitain seed-correctness, 
//once the sequence is stuck between the ambiguity from uniqu to repeatThreshold mode.
//Noted by KuanWeiLee 20180106
SeedFeature::SeedVector LongReadProbe::removeHitchhikingSeeds(SeedFeature::SeedVector initSeedVec)
{
	if(initSeedVec.size() < 2) return initSeedVec;
	
	for(SeedFeature::SeedVector::iterator iterQuery = initSeedVec.begin(); (iterQuery + 1) != initSeedVec.end(); iterQuery++)
	{
		SeedFeature& query = *iterQuery;
		SeedFeature::SeedVector::iterator iterSubject = iterQuery + 1;
		
		for(; iterSubject != initSeedVec.end(); iterSubject++)
		{
			SeedFeature& subject = *iterSubject;
			if((int)(subject.seedStartPos - query.seedEndPos) > m_params.radius) break;
			
			float freqDiff = (float)subject.maxFixedMerFreq/query.maxFixedMerFreq;
			
			subject.isHitchhiked |= (query.isRepeat && freqDiff < m_params.hhRatio);	//HIGH --> LOW
			query.isHitchhiked |= (subject.isRepeat && freqDiff > 1/m_params.hhRatio);	//LOW  --> HIGH
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
		std::ostream* pOutcastSeedWriter = createWriter(m_params.directory + "seed/error/" + readid + ".seed");
		*pOutcastSeedWriter << outcastSeedVec;
		delete pOutcastSeedWriter;
	}
	return finalSeedVec;
}
