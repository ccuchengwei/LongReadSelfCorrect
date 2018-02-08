/***************************/
/*** Seed Feature Body *****/
/***************************/
#include "SeedFeature.h"
#include "BWTAlgorithms.h"
#include "Util.h"
/*************************/
//Legacy code
SeedFeature::SeedFeature(
		size_t startPos,
		std::string str,
		bool repeat,
		size_t staticKmerSize,
		size_t repeatCutoff,
		size_t maxFixedMerFreq)
:	seedStr(str),
	seedStartPos(startPos),
	maxFixedMerFreq(maxFixedMerFreq),
	isRepeat(repeat),
	isHitchhiked(false),
	minKmerSize(staticKmerSize),
	freqUpperBound(repeatCutoff),
	freqLowerBound(repeatCutoff>>1)
{
	seedLength = seedStr.length();
	seedEndPos = seedStartPos + seedLength -1;
	startBestKmerSize = endBestKmerSize = staticKmerSize;
}
/*************************/
void SeedFeature::estimateBestKmerSize(const BWTIndexSet& indices)
{
	modifyKmerSize(indices, true);
	modifyKmerSize(indices, false);
}
//pole(true/false) ? start : end
//bit(1/-1) > 0 ? increase : decrease
void SeedFeature::modifyKmerSize(const BWTIndexSet& indices, bool pole)
{
	int& kmerSize = pole ? startBestKmerSize : endBestKmerSize;
	int& kmerFreq = pole ? startKmerFreq : endKmerFreq;
	const BWT* const pSelBWT = pole ? indices.pRBWT : indices.pBWT;
	std::string seed = pole ?  reverse(seedStr) : seedStr;
	kmerFreq = BWTAlgorithms::countSequenceOccurrences(seed.substr(seedLength - kmerSize), pSelBWT);
	int bit;
	if(kmerFreq > freqUpperBound)
		bit = 1;
	else if (kmerFreq < freqLowerBound)
		bit = -1;
	else
		return;
	const int freqBound = bit > 0 ? freqUpperBound : freqLowerBound;
	const int compFreqBound = bit > 0 ? freqLowerBound : freqUpperBound;
	const int sizeBound = bit > 0 ? sizeUpperBound : sizeLowerBound;
	
	while((bit^kmerFreq) > (bit^freqBound) && (bit^kmerSize) < (bit^sizeBound))
	{
		kmerSize += bit;
		kmerFreq = BWTAlgorithms::countSequenceOccurrences(seed.substr(seedLength - kmerSize), pSelBWT);
	}
	if((bit^kmerFreq) < (bit^compFreqBound))
	{
		kmerSize -= bit;
		kmerFreq = BWTAlgorithms::countSequenceOccurrences(seed.substr(seedLength - kmerSize), pSelBWT);
	}
}

