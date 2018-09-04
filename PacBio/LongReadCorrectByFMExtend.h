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
										size_t maxOverlap,
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
			inline bool isEmpty(){return m_leaves.empty();};

		// return size of leaves
			inline size_t size(){return m_leaves.size();};

		// return size of seed
			inline size_t getSeedSize(){return m_seedSize;};

		// return size of seed
			inline size_t getCurrentLength(){return m_currentLength;};

		size_t minTotalcount = 10000000;
		size_t totalcount = 0;

		size_t SelectFreqsOfrange(const size_t LowerBound, const size_t UpperBound, SONode3PtrList& newLeaves);

		std::pair<size_t,size_t> alnscore;
	private:

		//
		// Functions
		//
			void initialRootNode(const std::string& beginningkmer);
			void buildOverlapbyFMindex(IntervalTree<size_t>& fwdIntervalTree,IntervalTree<size_t>& rvcIntervalTree,const int& overlapSize);

			void extendLeaves(SONode3PtrList& newLeaves);
			void attempToExtend(SONode3PtrList& newLeaves, const bool isSuccessToReduce);
			void updateLeaves(SONode3PtrList& newLeaves,extArray& extensions,SAIOverlapNode3* leaf, const leafTable_t& lastPathInfo, leafTable_t& currPathInfo, size_t currLeavesNum);

			void refineSAInterval(SONode3PtrList& leaves, const size_t newKmerSize);

			int findTheBestPath(const SAIntervalNodeResultVector& results, FMWalkResult2& FMWResult);

			extArray getFMIndexExtensions(SAIOverlapNode3* leaf,const bool printDebugInfo);

			// prone the leaves without seeds in proximity
				bool PrunedBySeedSupport(SONode3PtrList& newLeaves);
			//Check if need reduce kmer size
				bool isInsufficientFreqs(SONode3PtrList& newLeaves);

			// Check if the leaves reach $
				bool isTerminated(SAIntervalNodeResultVector& results);

			bool isOverlapAcceptable(SAIOverlapNode3* currNode);
			bool isSupportedByNewSeed(SAIOverlapNode3* currNode, size_t smallSeedIdx, size_t largeSeedIdx);
			bool ismatchedbykmer(BWTInterval currFwdInterval,BWTInterval currRvcInterval);
			double computeErrorRate(SAIOverlapNode3* currNode);

			void printDebugData(debugPerExtArray& currDebugData);
			void printErrorRate(SONode3PtrList& currLeaves);

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
			const size_t m_maxOverlap;
			const BWT* m_pBWT;
			const BWT* m_pRBWT;
			const size_t m_PBcoverage;
			size_t m_min_SA_threshold;
			double m_errorRate;
			const size_t m_maxLeaves;
			const size_t m_seedSize;
			size_t m_repeatFreq;
			size_t m_localSimilarlykmerSize;
			const double m_PacBioErrorRate;

		// debug tools
			debugExtInfo m_Debug;

		size_t m_maxIndelSize;
		static thread_local std::vector<double> freqsOfKmerSize;

		// Optional parameters
			size_t m_maxfreqs;

		std::string m_query;
		size_t m_maxLength;
		size_t m_minLength;
		std::vector<BWTInterval> m_fwdTerminatedInterval;   //in rBWT
		std::vector<BWTInterval> m_rvcTerminatedInterval;   //in BWT

		SONode3PtrList m_leaves;
		SAIOverlapNode3* m_pRootNode;
		SONode3PtrList m_RootNodes;

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

		pathPack_t totalPathInfo;
};

#endif