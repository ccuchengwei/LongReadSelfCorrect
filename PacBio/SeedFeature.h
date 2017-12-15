#ifndef SeedFeature_H
#define SeedFeature_H

#include <iostream>
#include "BWTIndexSet.h"
class SeedFeature
{
	public:
		SeedFeature(size_t startPos, std::string str, bool repeat, size_t kmerSize, size_t repeatCutoff, size_t maxFixedMerFreqs=0);
		~SeedFeature(){};
		
		// append current seed string with extendedStr
		inline void append(std::string extendedStr, const SeedFeature& target)
		{
			seedStr += extendedStr;
			seedLength += extendedStr.length();
			//Upadate seed features in source to target
			startBestKmerSize = target.startBestKmerSize;
			endBestKmerSize = target.endBestKmerSize;
			isRepeat = target.isRepeat;
			maxFixedMerFreqs = target.maxFixedMerFreqs;
			seedStartPos = target.seedStartPos;
			seedEndPos = target.seedEndPos;
		};
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
		void estimateBestKmerSize(const BWTIndexSet& indices);
		
		size_t seedStartPos;
		size_t seedEndPos;
		size_t seedLength;
        size_t maxFixedMerFreqs;
		std::string seedStr;
		bool isRepeat;
		bool isHitchhiked;
		/*************************/
		//Unknown usage
		bool isPBSeed;
		bool isNextRepeat = false;
        bool isLargeVar = false;
		/*************************/
		// estimated by calling estimateBestKmerSize
		size_t startBestKmerSize;
		size_t endBestKmerSize;
		size_t startKmerFreq;
		size_t endKmerFreq;
		
	private:
		size_t minKmerSize;
		size_t freqUpperBound;
		size_t freqLowerBound;
		//size_t stepSize;
		void modifyKmerSize(const BWTIndexSet& indices, bool which);
		//estimate kmer size
		void increaseStartKmerSize(const BWT* pBWT);
		void decreaseStartKmerSize(const BWT* pBWT);

		//estimate kmer size
		void increaseEndKmerSize(const BWT* pBWT);
		void decreaseEndKmerSize(const BWT* pBWT);
};

#endif