//----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------
//
// Re-written from Jared Simpson's StringThreaderNode and StringThreader class
// The search tree represents a traversal through implicit FM-index graph
//
#ifndef OverlapTree_H
#define OverlapTree_H

#include <sstream>
#include "IntervalTree.h"
#include "SAINode.h"
#include "LongReadCorrtoolsByFMExtend.h"

class LongReadSelfCorrectByFMExtend
{
	public:
		LongReadSelfCorrectByFMExtend();
		LongReadSelfCorrectByFMExtend(
										const size_t caseNumber,
										const seedPair& extSeeds,
										const std::string& strBetweenSrcTarget,
										int m_disBetweenSrcTarget,
										size_t initkmersize,
										size_t maxResetSize,
										const FMextendParameters params,
										size_t m_min_SA_threshold = 3,
										const debugExtInfo debug = debugExtInfo(),
										double errorRate = 0.25,
										size_t repeatFreq = 256,
										size_t localSimilarlykmerSize = 100
									);

		~LongReadSelfCorrectByFMExtend();

		// extend all leaves one base pair during overlap computation
			int extendOverlap(FMWalkResult2& FMWResult);

		// return emptiness of leaves
			inline bool isEmpty(){return m_currLeavesKeptSize == 0;};

		// return size of leaves
			inline size_t size(){return m_currLeavesKeptSize;};

		// return size of match
			inline size_t getMatchedSize(){return m_matchedSize;};

		// return current length
			inline size_t getCurrentLength(){return m_currentLength;};

		size_t minTotalcount = 10000000;
		size_t totalcount = 0;

		size_t SelectFreqsOfrange(const size_t LowerBound, const size_t UpperBound, LeavesArray_t& newLeaves);

		std::pair<size_t,size_t> alnscore;
	private:

		//
		// Functions
		//
			void initialRootNode(const std::string& beginningkmer);
			void buildOverlapbyFMindex(IntervalTree<size_t>& fwdIntervalTree,IntervalTree<size_t>& rvcIntervalTree,const int& overlapSize);

			void extendLeaves();
			void attempToExtend(const debugBits_t& debugBits, std::string& debugStr, const leafTable_t& lastPathInfo, leafTable_t& currPathInfo);
			void updateLeaves(extArray& extensions,SAIOverlapNode3* leaf, const leafTable_t& lastPathInfo, leafTable_t& currPathInfo, size_t& currLeavesNum);

			void refineSAInterval(LeavesArray_t& leaves, const size_t newKmerSize);

			int findTheBestPath(const SAIntervalNodeResultVector& results, FMWalkResult2& FMWResult);

			extArray getFMIndexExtensions(SAIOverlapNode3* leaf, const debugBits_t& debugBits, std::stringstream& strBuffer);

			// prone the leaves without matches in proximity
				void PrunedByRelErrorRate(const leafTable_t& lastPathInfo);
				bool PrunedByAbsErrorRate();
			//Check if need reduce kmer size
				bool isInsufficientFreqs();

			// Check if the leaves reach $
				bool isTerminated(SAIntervalNodeResultVector& results);

			bool isOverlapAcceptable(SAIOverlapNode3* currNode);
			bool findNewMatch(SAIOverlapNode3* currNode, size_t smallMatchedIdx, size_t largeMatchedIdx);
			bool ismatchedbykmer(BWTInterval currFwdInterval,BWTInterval currRvcInterval);
			double computeErrorRate(SAIOverlapNode3* currNode);

			void printDebugData(debugPerExtArray& currDebugData, std::stringstream& strBuffer);
			void printErrorRate(LeavesArray_t& currLeaves);

		//
		// Data
		//

			const size_t m_case_number;
			size_t m_step_number;

			seedPair m_extSeeds;
			const std::string m_sourceSeed;
			const std::string m_strBetweenSrcTarget;
			const std::string m_targetSeed;
			const int m_disBetweenSrcTarget;
			const size_t m_initkmersize;
			const size_t m_minOverlap;
			const size_t m_maxResetSize;
			const BWT* m_pBWT;
			const BWT* m_pRBWT;
			const size_t m_PBcoverage;
			size_t m_min_SA_threshold;
			double m_errorRate;
			const size_t m_maxLeaves;
			const size_t m_matchedSize;
			size_t m_repeatFreq;
			size_t m_localSimilarlykmerSize;
			const double m_PacBioErrorRate;
			double m_minMatchedErrorRate;

		// debug tools
			debugExtInfo m_Debug;
			short int    fieldNum;

		size_t m_maxIndelSize;
		static thread_local std::vector<double> freqsOfKmerSize;

		// Optional parameters
			size_t m_maxfreqs;

		std::string m_query;
		size_t m_maxLength;
		size_t m_minLength;
		std::vector<BWTInterval> m_fwdTerminatedInterval;   //in rBWT
		std::vector<BWTInterval> m_rvcTerminatedInterval;   //in BWT

		// Node Info
			// Leaves
				LeavesArray_t m_currLeaves;
				LeavesArray_t m_nextLeaves;
			// Leaves Number
				// Current leaves containing pruned leaves.
					size_t m_currLeavesTotalSize;
				// Current leaves not containing pruned leaves.
					size_t m_currLeavesKeptSize;
				// The next leaves containing pruned leaves.
					size_t m_nextLeavesTotalSize;
				// The next leaves not containing pruned leaves.
					size_t m_nextLeavesKeptSize;
			// Roots
				SAIOverlapNode3* m_pRootNode;
				LeavesArray_t m_RootNodes;
			// path
				pathPack_t totalPathInfo;
				const bool isSeparatedInitKmer;
				size_t m_pathOffset;

		size_t m_currentLength;
		size_t m_currentKmerSize;

		IntervalTree<size_t> m_fwdIntervalTree;
		IntervalTree<size_t> m_rvcIntervalTree;

		IntervalTree<size_t> m_fwdIntervalTree2;
		IntervalTree<size_t> m_rvcIntervalTree2;

		size_t RemainedMaxLength;
		size_t m_totalExtNum;
		size_t m_currTotalExtNum;
		size_t m_currExtNum;
		size_t m_accumLeaves;

};

#endif