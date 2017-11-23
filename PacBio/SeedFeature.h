#ifndef SeedFeature_H
#define SeedFeature_H

#include <iostream>
#include "BWTIndexSet.h"
struct SeedFeature
{
	public:
		SeedFeature(size_t startPos, std::string str, bool repeat, size_t kmerSize, size_t repeatCutoff);

		SeedFeature(){};
		~SeedFeature(){};
		
		// append current seed string with extendedStr
		inline void append(std::string extendedStr)
		{
			seedStr += extendedStr;
			seedLength += extendedStr.length();
			seedStartPos += extendedStr.length();
			seedEndPos += extendedStr.length();
		};
		inline bool isSmall(){return seedLength<=17?true:false ;};
		inline void setBestKmerSize(size_t staticKmerSize)
		{
			startBestKmerSize = endBestKmerSize = staticKmerSize;;
		};
		void estimateBestKmerSize(const BWT* pBWT);
		
		size_t seedStartPos;
		size_t seedEndPos;
		size_t seedLength;
        size_t maxFixedMerFreqs;
		std::string seedStr;
		bool isRepeat;
		bool isPBSeed;
		bool isNextRepeat = false;
        bool isLargeVar = false;
		
		// estimated by calling estimateBestKmerSize
		size_t startBestKmerSize;
		size_t endBestKmerSize;
		size_t startKmerFreq;
		size_t endKmerFreq;
		
	private:
		size_t freqUpperBound;
		size_t freqLowerBound;
		size_t minKmerSize;
		size_t stepSize;
		//estimate kmer size
		void increaseStartKmerSize(const BWT* pBWT);
		void decreaseStartKmerSize(const BWT* pBWT);

		//estimate kmer size
		void increaseEndKmerSize(const BWT* pBWT);
		void decreaseEndKmerSize(const BWT* pBWT);
};
#endif