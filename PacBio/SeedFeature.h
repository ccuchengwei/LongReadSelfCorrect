#ifndef SeedFeature_H
#define SeedFeature_H

#include <iostream>
#include "BWTIndexSet.h"
class SeedFeature
{
	public:
		/*************************/
		//Legacy code
		SeedFeature(size_t startPos, std::string str, bool repeat, size_t kmerSize, size_t repeatCutoff, size_t maxFixedMerFreq=0);
		/*************************/
		SeedFeature(
				std::string str,
				size_t startPos,
				size_t frequency,
				bool repeat,
				size_t kmerSize,
				size_t PBcoverage)
		:	seedStr(str),
			seedLength(seedStr.length()),
			seedStartPos(startPos),
			seedEndPos(startPos + seedLength - 1),
			maxFixedMerFreq(frequency),
			isRepeat(repeat),
			isHitchhiked(false),
			startBestKmerSize(kmerSize),
			endBestKmerSize(kmerSize),
			sizeUpperBound(seedLength),
			sizeLowerBound(kmerSize),
			freqUpperBound(PBcoverage >> 1),
			freqLowerBound(PBcoverage >> 2){ }
				
		~SeedFeature(){ }
		
		// append current seed string with extendedStr
		inline void append(std::string extendedStr, const SeedFeature& target)
		{
			seedStr += extendedStr;
			seedLength += extendedStr.length();
			//Upadate seed features in source to target
			startBestKmerSize = target.startBestKmerSize;
			endBestKmerSize = target.endBestKmerSize;
			isRepeat = target.isRepeat;
			maxFixedMerFreq = target.maxFixedMerFreq;
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
		
		std::string seedStr;
		size_t seedLength;
		size_t seedStartPos;
		size_t seedEndPos;
        size_t maxFixedMerFreq;
		bool isRepeat;
		bool isHitchhiked;
		/*************************/
		//Unknown usage
		bool isPBSeed;
		bool isNextRepeat = false;
        bool isLargeVar = false;
		/*************************/
		// estimated by calling estimateBestKmerSize
		int startBestKmerSize;
		int endBestKmerSize;
		int startKmerFreq;
		int endKmerFreq;
		
	private:
		size_t minKmerSize;
		int sizeUpperBound;
		int sizeLowerBound;
		int freqUpperBound;
		int freqLowerBound;
		//size_t stepSize;
		void modifyKmerSize(const BWTIndexSet& indices, bool which);
		/*
		//estimate kmer size
		void increaseStartKmerSize(const BWT* pBWT);
		void decreaseStartKmerSize(const BWT* pBWT);

		//estimate kmer size
		void increaseEndKmerSize(const BWT* pBWT);
		void decreaseEndKmerSize(const BWT* pBWT);
		*/
};

#endif