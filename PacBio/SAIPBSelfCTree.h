//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------

#ifndef SAIPBSelfCTree_H
#define SAIPBSelfCTree_H

#include <list>

#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SAINode.h"
#include "HashMap.h"




struct END
{
	std::string endingKmer;
	int maxLength;
	int minLength;
    BWTInterval fwdTerminatedInterval;   //in rBWT
    BWTInterval rvcTerminatedInterval;   //in BWT
};

// class for storing various features of each kmer in the hashtable, including pos, freq, ...
class KmerFeatures{
	private:
		std::vector<long long int> m_sumOfFreq;
		std::vector<long long int> m_sumOfPos;
		long long int m_intervalSize;
		long long int m_totalFreq;
		long long int m_totalSum;
		double m_maxAvgFreq;
		bool m_isVisited;

    public: 
		// maxIntervalSize is overestimated to be the max length
		KmerFeatures(long long int pos, size_t maxIntervalSize, size_t intervalSize=35):
			m_intervalSize(intervalSize), 
			m_totalFreq(0),m_totalSum(0), m_maxAvgFreq(0), m_isVisited(false)
			{
				m_sumOfFreq.resize(maxIntervalSize/intervalSize + 1);
				m_sumOfPos.resize(maxIntervalSize/intervalSize + 1);
				add(pos);
			}
			
		~KmerFeatures(){
		}

		// increment frequence and compute approximate pos
		void add(long long int pos)
		{
			m_totalFreq++;
			m_totalSum+=pos;
			
			// compute sampled index
			int index = pos/m_intervalSize;
			// extension may exceed the expected length
           
			if(index < 0)
				index = 0;
			else if(index > (int) m_sumOfFreq.size()-1)
				index = m_sumOfFreq.size()-1;

			
			m_sumOfFreq[index]++;
			m_sumOfPos[index]+=pos;
		}
		
		long long int getTotalFreq()
		{
			return m_totalFreq;
		}

		long long int getTotalSum()
		{
			return m_totalSum;
		}

		// convert pos into index in the m_sumOfPos
		long long int getSumOfFreq(long long int pos)
		{
			int index = pos/m_intervalSize;
           
			long long int sumOfFreq = m_sumOfFreq[index];
			if(index > 0)
				sumOfFreq += m_sumOfFreq[index-1];
			if(index < (int) m_sumOfFreq.size()-1)
				sumOfFreq += m_sumOfFreq[index+1];
				
			return sumOfFreq;
		}
		
		// convert pos into index in the m_sumOfPos
		long long int getSumOfPos(long long int pos)
		{
			int index = pos/m_intervalSize;
			long long int sumOfPos = m_sumOfPos[index];
			if(index > 0)
				sumOfPos += m_sumOfPos[index-1];
			if(index < (int) m_sumOfPos.size()-1)
				sumOfPos += m_sumOfPos[index+1];
				
			return sumOfPos;
		}
		
		void printFreq()
		{
			for(size_t i=0; i<m_sumOfFreq.size(); i++)
			{
				std::cout << i <<":" << m_sumOfFreq[i] << "\t" << m_sumOfPos[i] <<"\n";
			}
		}
		
		inline void setVisited(bool visited)
		{
			m_isVisited = visited;
		}
		
		inline bool isVisited()
		{
			return m_isVisited;
		}
		
		inline void setMaxAvgFreq(double newAvgFreq)
		{
			m_maxAvgFreq = newAvgFreq;
		}
		
		inline double getMaxAvgFreq()
		{
			return m_maxAvgFreq;
		}
};

class SAIPBSelfCorrectTree
{
    public: 
		SAIPBSelfCorrectTree(const BWT* pBWT, const BWT* pRBWT,const std::string rawSeq,const size_t srcmaxLength, size_t min_SA_threshold=2, int m_maxLeavesAllowed=64);
		
        ~SAIPBSelfCorrectTree();

		// Collect kmers from overlapping reads via right extension of single seed
		// void addHashFromSingleSeedUsingFMExtension(std::string& seedStr, size_t hashKmerSize, size_t maxLength);
		
		// LF-mapping of each SA index in the interval independently using loop instead of BFS tree expansion
		// Contaminated reads often are simple repeats C* or T* with large freq
		// Return false if the freq is way too large
		size_t addHashBySingleSeed(std::string& seedStr, size_t largeKmerSize, size_t smallKmerSize, size_t maxLength, bool skipRepeat, int expectedLength=-1);

		size_t addHashByExtSeed(std::string& extStr, size_t srcLength, size_t largeKmerSize, size_t smallKmerSize, size_t extLength, size_t maxLength);

		// Collect kmers from overlapping reads via FM-index walk between two seeds
		// int addHashFromPairedSeed(std::string seedStr, std::vector<std::pair<int, std::string> >  &targets, size_t hashKmerSize);
		
        // find feasible extension using kmers collected by addHashFromSingleSeed or addHashFromPairedSeed
		int mergeTwoSeedsUsingHash(const std::string &src, const std::string &dest, std::string &mergedseq, 
									size_t minOverlap, size_t maxLeaves, size_t minLength, size_t maxLength, size_t expectedLength);
		
		int extOverUnique(const std::string &src, const std::string &dest, std::string &mergedseq, 
									size_t minOverlap, size_t maxLeaves, size_t minLength, size_t maxLength, size_t expectedLength);

		int extOverRepeat(const std::string &src, const std::string &dest, std::string &mergedseq, 
									size_t minOverlap, size_t maxLeaves, size_t minLength, size_t maxLength, size_t expectedLength,std::string currseq);
									
		void collectExtSeeds();
		
        // Print all the strings represented by the tree
        void printAll();
		void printLeaves(size_t hashKmerSize);
        size_t hashkmerfreqs(std::string fwdkmer,size_t kemrposition);
    private:

        //
        // Functions
        //
		// create root node using src
		void initializeSearchTree(std::string src, size_t hashKmerSize);
		void initializeTerminalIntervals(std::string dest, size_t hashKmerSize);

		void attempToExtendUsingHash(STNodePtrList &newLeaves, size_t hashKmerSize, size_t minExtFreq);
        void refineSAInterval(size_t newKmer);
		std::vector<std::pair<std::string, BWTIntervalPair> > getFMIndexRightExtensions(SAIntervalNode* pNode, const size_t IntervalSizeCutoff);
		
		void insertKmerToHash(std::string& insertedKmer, size_t seedStrLen, size_t currentLength, size_t smallKmerSize, size_t maxLength, int expectedLength);

		bool isExtensionValid(std::string fwdkmer, double& currAvgFreq, size_t& kmerFreq,size_t bcount);
		
		void NewextendByNewSeeds(const std::string& src, const std::string& dest, size_t hashKmerSize, size_t maxLength,std::string currseq);
		
		void extendByNewSeeds(const std::string& src, const std::string& dest, size_t hashKmerSize, size_t maxLength);
		
		void rollBackExtension(size_t newKmerSize, bool btrimAll = false);
		
        // Check if the leaves can be extended no further
		bool isTerminated(SAIntervalNodeResultVector& results, END target);
		bool isTerminated(SAIntervalNodeResultVector& results);

		bool isLowComplexity (std::string seq, float threshold = 0.8, bool checkTwoBase = false)
		{
			size_t seqLen = seq.length();
			size_t countG =0 ;
			size_t countC =0 ;
			size_t countT =0 ;
			size_t countA =0 ;

			size_t dimerCount = 0;
			char prevChar = ' ';
			for (size_t i=0; i<seqLen; i++)
			{
				switch(seq[i]){
					case 'A': countA ++ ;break;
					case 'T': countT ++ ;break;
					case 'C': countC ++ ;break;
					case 'G': countG ++ ;break;
					default:  assert(false);
				}
				
				if(seq[i] == prevChar)
					dimerCount++;
					
				prevChar = seq[i];
			}

			bool isOneBaseLowCom = ((float) countA/seqLen >= threshold ) || ((float) countT/seqLen >= threshold)
					|| ((float) countC/seqLen >= threshold ) || ((float) countG/seqLen >= threshold);
					
			// float countAC = countA+countC;
			// float countAG = countA+countG;
			// float countAT = countA+countT;
			// float countCG = countC+countG;
			// float countGT = countG+countT;
			// float countCT = countC+countT;
			
			// bool isTwoBaseLowCom = countAC/seqLen >= threshold || countAG/seqLen >= threshold ||
									// countAT/seqLen >= threshold || countCG/seqLen >= threshold ||
									// countGT/seqLen >= threshold || countCT/seqLen >= threshold;
									
			if ( isOneBaseLowCom || (checkTwoBase && (float)dimerCount/(seqLen-1)>0.68) )
			{
				// std::cout << countA << "\t" << dimerCount << "\n";
				return true;
			}
			return false;
		}
        //
        // Data
        //
        
		int findBestExtSeq(SAIntervalNodeResultVector& results, const std::string &src, const std::string &dest, std::string &mergedseq,size_t hashKmerSize, size_t expectedLength, size_t maxLength, size_t maxLeaves);
		
        const BWT* m_pBWT;
        const BWT* m_pRBWT;
        const std::string m_rawSeq;
        const size_t m_maxLength;
        size_t m_min_SA_threshold;
		size_t m_maxLeavesAllowed;
		
		int m_expectedLength;
		int m_currentLength;
		int m_seedLength;
		
		std::vector<END> ENDs;

        SAIntervalNode* m_pRootNode;
        STNodePtrList m_leaves;

        BWTInterval m_fwdTerminatedInterval;   //in rBWT
        BWTInterval m_rvcTerminatedInterval;   //in BWT
		
		// multiple means for repeats leading to multinomial occurrences
		DenseHashMap<std::string, KmerFeatures*, StringHasher> kmerHash;
		typedef DenseHashMap<std::string, KmerFeatures*> :: iterator kmerHashiter;
		
		bool m_isSourceRepeat;
		bool m_isTargetRepeat;
		size_t m_sourceKmerSize;
		size_t m_targetKmerSize;
		
		bool m_isLargeLeaveRemoved;
		size_t m_prevExtPos;
		// bool debug;
};

#endif
