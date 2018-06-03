#ifndef SeedFeature_H
#define SeedFeature_H

#include <iostream>
#include <map>
#include "BWTIndexSet.h"
class SeedFeature
{
	public:
		
		typedef std::vector<SeedFeature> SeedVector;
		static std::map<std::string, SeedVector>& Log();
		friend std::ostream& operator<<(std::ostream& out, const SeedVector& vec);
		
		SeedFeature(std::string str, int startPos, int frequency, bool repeat, int kmerSize, int PBcoverage);
		~SeedFeature(void) = default;
		
		// adjust start/end kmer for future FMWalk
		void estimateBestKmerSize(const BWTIndexSet& indices);
		
		// append current seed string with extendedStr
		inline void append(std::string extendedStr, const SeedFeature& target)
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
		
		std::string seedStr;
		int seedLen;
		int seedStartPos;
		int seedEndPos;
        int maxFixedMerFreq;
		bool isRepeat;
		bool isHitchhiked;
		int startBestKmerSize;
		int endBestKmerSize;
		int startKmerFreq;
		int endKmerFreq;
		
		//Legacy part
		/***********/
		SeedFeature(size_t startPos, std::string str, bool repeat, size_t kmerSize, size_t repeatCutoff, size_t maxFixedMerFreq=0);
		inline void append(std::string extendedStr)
		{
			seedStr += extendedStr;
			seedLen += extendedStr.length();
			seedStartPos += extendedStr.length();
			seedEndPos += extendedStr.length();
		};
		inline bool isSmall(){return seedLen<=17?true:false ;};
		inline void setBestKmerSize(size_t staticKmerSize)
		{
			startBestKmerSize = endBestKmerSize = staticKmerSize;;
		}
		
		bool isPBSeed;
		bool isNextRepeat = false;
        bool isLargeVar = false;
		int minKmerSize;
		/***********/
	private:
		int sizeUpperBound;
		int sizeLowerBound;
		int freqUpperBound;
		int freqLowerBound;
		void modifyKmerSize(const BWTIndexSet& indices, bool which);
};

#endif