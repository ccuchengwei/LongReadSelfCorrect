/***************************/
/*** Seed Feature Body *****/
/***************************/
#include "SeedFeature.h"
#include "BWTAlgorithms.h"
#include "Util.h"
SeedFeature::SeedFeature(size_t startPos, std::string str, bool repeat, size_t staticKmerSize, size_t repeatCutoff, size_t maxFixedMerFreqs)
	:seedStartPos(startPos), maxFixedMerFreqs(maxFixedMerFreqs), seedStr(str), isRepeat(repeat), isHitchhiked(false), 
	minKmerSize(staticKmerSize), freqUpperBound(repeatCutoff),freqLowerBound(repeatCutoff>>1)
{
	seedLength = seedStr.length();
	seedEndPos = seedStartPos + seedLength -1;
	startBestKmerSize = endBestKmerSize = staticKmerSize;
}

void SeedFeature::estimateBestKmerSize(const BWTIndexSet& indices)
{
	modifyKmerSize(indices, true);
	modifyKmerSize(indices, false);
}
//which(true/false) ? start : end
//type(1/-1) > 0 ? increase : decrease
void SeedFeature::modifyKmerSize(const BWTIndexSet& indices, bool which)
{
	size_t& kmerFreq = which ? startKmerFreq : endKmerFreq;
	size_t& kmerSize = which ? startBestKmerSize : endBestKmerSize;
	const BWT* pSelBWT = which ? indices.pRBWT : indices.pBWT;
	std::string seed = which ?  reverse(seedStr) : seedStr;
	kmerFreq = BWTAlgorithms::countSequenceOccurrences(seed.substr(seedLength - kmerSize), pSelBWT);
	int type;
	if(kmerFreq > freqUpperBound)
		type = 1;
	else if (kmerFreq < freqLowerBound)
		type = -1;
	else
		return;
	size_t freqBound = type > 0 ? freqUpperBound : freqLowerBound;
	size_t compFreqBound = type > 0 ? freqLowerBound : freqUpperBound;
	size_t sizeBound = type > 0 ? seedLength : minKmerSize;
	
	while((type^kmerFreq) > (type^freqBound) && (type^kmerSize) < (type^sizeBound))
	{
		kmerSize += type;
		kmerFreq = BWTAlgorithms::countSequenceOccurrences(seed.substr(seedLength - kmerSize), pSelBWT);
	}
	if(type*kmerFreq < type*compFreqBound)
	{
		kmerSize -= type;
		kmerFreq = BWTAlgorithms::countSequenceOccurrences(seed.substr(seedLength - kmerSize), pSelBWT);
	}
}

