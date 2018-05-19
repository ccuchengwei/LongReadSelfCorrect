#include "SeedFeature.h"
#include "BWTAlgorithms.h"
#include "Util.h"

//Legacy part
/***********/
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
	seedLen = seedStr.length();
	seedEndPos = seedStartPos + seedLen -1;
	startBestKmerSize = endBestKmerSize = staticKmerSize;
}
/***********/

std::map<std::string, SeedFeature::SeedVector>& SeedFeature::Log()
{
	static std::map<std::string, SeedVector> log;
	return log;
}

void SeedFeature::write(std::ostream& out, const SeedVector& vec)
{
	for(auto& iter : vec)
		out
		<< iter.seedStr << '\t'
		<< iter.maxFixedMerFreq << '\t' 
		<< iter.seedStartPos << '\t'
		<< (iter.isRepeat ? "Yes" : "No") << '\n';
}

SeedFeature::SeedFeature(
		std::string str,
		int startPos,
		int frequency,
		bool repeat,
		int kmerSize,
		int PBcoverage)
:	seedStr(str),
	seedLen(seedStr.length()),
	seedStartPos(startPos),
	seedEndPos(startPos + seedLen - 1),
	maxFixedMerFreq(frequency),
	isRepeat(repeat),
	isHitchhiked(false),
	startBestKmerSize(kmerSize),
	endBestKmerSize(kmerSize),
	sizeUpperBound(seedLen),
	sizeLowerBound(kmerSize),
	freqUpperBound(PBcoverage >> 1),
	freqLowerBound(PBcoverage >> 2){ }

void SeedFeature::append(std::string extendedStr, const SeedFeature& target)
{
	seedStr += extendedStr;
	seedLen += extendedStr.length();
	//Upadate seed features of source into target
	startBestKmerSize = target.startBestKmerSize;
	endBestKmerSize = target.endBestKmerSize;
	isRepeat = target.isRepeat;
	maxFixedMerFreq = target.maxFixedMerFreq;
	seedStartPos = target.seedStartPos;
	seedEndPos = target.seedEndPos;
}
		
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
	kmerFreq = BWTAlgorithms::countSequenceOccurrences(seed.substr(seedLen - kmerSize), pSelBWT);
	int bit;
	if(kmerFreq > freqUpperBound)
		bit = 1;
	else if (kmerFreq < freqLowerBound)
		bit = -1;
	else
		return;
	const int freqBound     = bit > 0 ? freqUpperBound : freqLowerBound;
	const int corsFreqBound = bit > 0 ? freqLowerBound : freqUpperBound;
	const int sizeBound = bit > 0 ? sizeUpperBound : sizeLowerBound;
	
	while((bit^kmerFreq) > (bit^freqBound) && (bit^kmerSize) < (bit^sizeBound))
	{
		kmerSize += bit;
		kmerFreq = BWTAlgorithms::countSequenceOccurrences(seed.substr(seedLen - kmerSize), pSelBWT);
	}
	if((bit^kmerFreq) < (bit^corsFreqBound))
	{
		kmerSize -= bit;
		kmerFreq = BWTAlgorithms::countSequenceOccurrences(seed.substr(seedLen - kmerSize), pSelBWT);
	}
}

