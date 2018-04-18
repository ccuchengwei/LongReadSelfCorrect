//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------

//
// Re-written from Jared Simpson's StringThreaderNode and StringThreader class
// The search tree represents a traversal through implicit FM-index graph
//
#ifndef SAINTERVALTREE_H
#define SAINTERVALTREE_H

#include <list>
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SAINode.h"

class LongReadExtend
{
	public:
		LongReadExtend(
				const std::string& _source,
				const std::string& _target,
				int _minOverlap,
				int _maxOverlap,
				int _maxLength,
				int _maxLeaves,
				BWTIndexSet _indices,
				int _SA_threshold = 3);
		
		~LongReadExtend();

		//return the merged string
		int bindTwoSeeds(std::string &mergedseq);
		
		/*
		// validate if each kmer in the read is suuported by at least minOverlap overlap 
		size_t getKmerCoverage(){return maxKmerCoverage;};
		size_t getMaxUsedLeaves(){return maxUsedLeaves;};
		bool isBubbleCollapsed(){return isBubbleCollapsed;}

		// Print all the strings represented by the tree
		void printAll();
		*/
	private:

		//
		// Functions
		//
		void extendLeaves();
		void attempToExtend(SAINodePtrList &newLeaves);
		std::vector<std::pair<std::string, BiBWTInterval> > getFMIndexExtensions(SAIPBNode* pNode);
		void refineSAInterval(int newKmerSize);

		// Check if the leaves can be extended no further
	//	bool isTerminated(SAIntervalNodeResultVector& results);
		bool isTerminated(std::vector<std::string>& results);
	//	bool isTwoReadsOverlap(std::string & mergedseq);
		/*
		size_t calculateKmerCoverage (const std::string & seq , size_t kmerLength , const BWT* pBWT);
		bool replaceLowFreqKmer (std::string & seq , size_t kmerLength);

		void removeLeavesByRepeatKmer();
		*/
		//
		// Data
		//
		const std::string source;
		const std::string target;
		int minOverlap;
		int maxOverlap;
		int maxLength;
		int maxLeaves;
		BWTIndexSet indices;
		int min_SA_threshold;

		SAIPBNode* pRootNode;
		SAINodePtrList leaves;

		int currentLength;
		int currentKmerSize;
	//	size_t maxKmerCoverage;
	//	size_t maxUsedLeaves;
		bool isBubbleCollapsed;
		
		BiBWTInterval destInterval;
};

#endif
