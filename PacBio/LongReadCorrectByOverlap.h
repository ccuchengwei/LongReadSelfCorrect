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

#include <list>
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SAINode.h"
#include "IntervalTree.h"
#include "HashtableSearch.h"
#include "SeedFeature.h"

// Define int type.
	typedef int64_t kmerFreq_t;

struct FMWalkResult2
{
	std::string mergedSeq;
	int alnScore;
	double kmerFreq;
};

struct FMextendParameters
{
	public:
		FMextendParameters(BWTIndexSet indices, int idmerLength, int maxLeaves,int minKmerLength,size_t PBcoverage, double ErrorRate):
			indices(indices),
			idmerLength(idmerLength),
			maxLeaves(maxLeaves),
			minKmerLength(minKmerLength),
			PBcoverage(PBcoverage),
			ErrorRate(ErrorRate){};

		FMextendParameters(){};
		BWTIndexSet indices;
		int idmerLength;
		int maxLeaves;
		int minKmerLength;
		size_t PBcoverage;
		double ErrorRate;

};

struct simpleSeedFeature
{
	public:
		simpleSeedFeature(const SeedFeature& seed, bool isHead):
			seq     (seed.seedStr),
			end     (seed.seedEndPos),
			isRepeat(seed.isRepeat)
			{
				if (isHead)
					start = seed.seedEndPos - seed.seedLen+1;
				else
					start = seed.seedStartPos;
			};

		std::string seq;
		size_t start;
		size_t end;
		bool   isRepeat;
};

struct seedPair
{
	public:
		seedPair(
					const SeedFeature& source,
					const SeedFeature& target,
					bool isPosStrand = true
				):
					source(source,true),
					target(target,false),
					isPosStrand(isPosStrand) {};

		void reverseSeed()
			{
				std::swap(source,target);
				source.seq = reverseComplement(source.seq);
				target.seq = reverseComplement(target.seq);
				isPosStrand = !isPosStrand;
			}

		void reduceSourceBy(size_t length)
			{
				source.seq = source.seq.substr(length);
				if (isPosStrand)
				{
					source.start = source.start + length;
				}
				else
				{
					source.end   = source.end   - length;
				}
			}

		simpleSeedFeature  source;
		simpleSeedFeature  target;
		bool    isPosStrand;
};

struct forceExtInfo
{
	forceExtInfo():
					isSpecifiedPath(false),
					specifiedPath(""){};

	void setSpecPath(
						const bool& isSpecifiedPath,
						const std::string& specifiedPath
					)
					{
						this -> isSpecifiedPath = isSpecifiedPath;
						this -> specifiedPath   = specifiedPath;
					};

	bool isEmpty()
	{
		return specifiedPath.empty();
	}

	bool isMatchBase(const std::string& currBase, const size_t position)
	{
		return (currBase.back() == specifiedPath.at(position));
	}

	bool isMatchBase(const char currBase, const size_t position)
	{
		return (currBase == specifiedPath.at(position));
	}

	bool isIllegalLen(size_t position)
	{
		return (position >= specifiedPath.length());
	}

	void RCSeq()
	{
		if (!isEmpty())
			specifiedPath = reverseComplement(specifiedPath);
	}

	std::string getSpecifiedPath()
	{
		return specifiedPath;
	}

		bool isSpecifiedPath;

	private:
		std::string specifiedPath;

};

struct debugExtInfo
{
	public:
		debugExtInfo(
						bool isDebug = false,
						std::ostream* debug_file     = nullptr,
						std::ostream* debug_finalSeq = nullptr,
						std::string readID = "",
						int  caseNum = 0,
						forceExtInfo ref = forceExtInfo()
					):
						isDebug(isDebug),
						debug_file(debug_file),
						debug_finalSeq(debug_finalSeq),
						readID(readID),
						caseNum(caseNum),
						ref(ref) {};

		const bool isDebug;
		std::ostream* debug_file;
		std::ostream* debug_finalSeq;
		std::string   readID;
		size_t        caseNum;
		forceExtInfo  ref;

};

struct FMidx_t
{
	public:
		FMidx_t  (
					const std::string& s,
					const BWTInterval& fwdInterval,
					const BWTInterval& rvcInterval
				):
					SearchLetters(s),
					fwdInterval(fwdInterval),
					rvcInterval(rvcInterval),
					kmerFrequency(fwdInterval.size() + rvcInterval.size())
				{};
		FMidx_t   (
					const char c,
					const BWTInterval& fwdInterval,
					const BWTInterval& rvcInterval
				):
					SearchLetters(std::string(1,c)),
					fwdInterval(fwdInterval),
					rvcInterval(rvcInterval),
					kmerFrequency(fwdInterval.size() + rvcInterval.size())
				{};
		void setInterval(const BWTInterval& fwdInterval,const BWTInterval& rvcInterval)
			{
				this -> fwdInterval   = fwdInterval;
				this -> rvcInterval   = rvcInterval;
				this -> kmerFrequency = fwdInterval.size() + rvcInterval.size();
			}

		BWTInterval getFwdInterval()
			{
				return fwdInterval;
			}

		BWTInterval getRvcInterval()
			{
				return rvcInterval;
			}

		kmerFreq_t getKmerFrequency()
			{
				return kmerFrequency;
			}

		char getLetter()
			{
				return SearchLetters[0];
			}

		bool isVaildIntervals()
			{
				return (fwdInterval.isValid() || rvcInterval.isValid());
			}

		const std::string SearchLetters;

	private:
		BWTInterval fwdInterval;
		BWTInterval rvcInterval;
		kmerFreq_t kmerFrequency;

};
typedef std::vector<FMidx_t> extArray;

struct leafInfo
{
	public:
		leafInfo(SAIOverlapNode3* leafNode, const size_t lastLeafNum):leafNodePtr(leafNode), lastLeafID(lastLeafNum)
			{
				const std::string leafLabel = leafNode -> getFullString();
				tailLetterCount   = 0;

				for(auto reverseIdx = leafLabel.crbegin(); reverseIdx != leafLabel.crend(); ++reverseIdx)
				{
					std::string suffixLetter(1,(*reverseIdx));
					if (reverseIdx == leafLabel.crbegin())
						tailLetter = suffixLetter;
					if (tailLetter == suffixLetter)
						tailLetterCount++;
					else
						break;
				}

				kmerFrequency =   (leafNode -> fwdInterval).size()
								+ (leafNode -> rvcInterval).size();
			}
		leafInfo(SAIOverlapNode3* currNode, const leafInfo& leaf, FMidx_t& extension, const size_t currLeavesNum)
			{
				const std::string& extLabel = extension.SearchLetters;

				// Set kmerFrequency
					kmerFrequency = extension.getKmerFrequency();

				// Set currNode
					// Copy the intervals
						currNode->fwdInterval = extension.getFwdInterval();
						currNode->rvcInterval = extension.getRvcInterval();
						currNode->addKmerCount( kmerFrequency );
					// currOverlapLen/queryOverlapLen always increase wrt each extension
					// in order to know the approximate real-time matched length for terminal/containment processing
						currNode -> currOverlapLen++;
						currNode -> queryOverlapLen++;
					leafNodePtr = currNode;

				// Set lastLeafID
					lastLeafID = currLeavesNum;

				// Set tailLetter and its counter
					if (leaf.tailLetter == extLabel)
					{
						tailLetter      = leaf.tailLetter;
						tailLetterCount = leaf.tailLetterCount + 1;
					}
					else
					{
						tailLetter = extLabel;
						tailLetterCount = 1;
					}
			}

		SAIOverlapNode3* leafNodePtr;
		size_t lastLeafID;
		kmerFreq_t kmerFrequency;

		std::string tailLetter;
		size_t tailLetterCount;

};
typedef std::list<leafInfo> leafList;

struct dominantBase
{
	public:
		dominantBase(const std::string& goalBase):
			identicalBase(goalBase[0]) {};

		int getMaxFreq()
			{
				return freqs[0];
			}

		void setFreq(FMidx_t currExt)
			{
				int  currFreq = currExt.getKmerFrequency();
				char currBase = currExt.getLetter();

				if( currFreq >  freqs[0]
				|| (currFreq == freqs[0] && currBase == identicalBase))
				{
					std::swap(freqs[0], currFreq);
					std::swap(bases[0], currBase);
				}
				if( currFreq >  freqs[1]
				|| (currFreq == freqs[1] && currBase == identicalBase))
				{
					std::swap(freqs[1], currFreq);
					std::swap(bases[1], currBase);
				}
			}

		bool isDominant()
			{
				return (bases[0] == identicalBase || bases[1] == identicalBase);
			}

		int  freqs[2] = {0,0};
		char bases[2] = {0,0};
		char identicalBase = 0;
};

struct match_t
{
	match_t   (
					bool isMatched=false,
					kmerFreq_t totalNum  = 0,
					kmerFreq_t matchFreq = 0
				):
					isMatched(isMatched),
					totalNum(totalNum),
					matchFreq(matchFreq) {};

		const bool isMatched;
		const kmerFreq_t totalNum;
		const kmerFreq_t matchFreq;
};
typedef std::vector<match_t> MatchArray;

struct debugPerExtInfo
{
	debugPerExtInfo (
						FMidx_t currExt,
						match_t matchInfo = 0
					):
						matchInfo(matchInfo)
					{
						(this  -> kmerFreq ) = currExt.getKmerFrequency();
						(this -> errorType ).reserve(4);
						(this  -> kmerRatio) = 0;
					}

	void setErrorData   (
							const double       kmerRatio,
							const std::string& errorString
						)
					{
						(this -> kmerRatio) = kmerRatio;
						(this -> errorType) = "X" + errorString;
					}
	void extSuccessSet()
					{
						(this -> errorType).front() = 'O';
					}

	std::string errorType;
	kmerFreq_t  kmerFreq;
	double      kmerRatio;
	match_t     matchInfo;
};
typedef std::vector<debugPerExtInfo> debugPerExtArray;

class LongReadSelfCorrectByOverlap
{
	public:
		LongReadSelfCorrectByOverlap();
		LongReadSelfCorrectByOverlap(
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

        ~LongReadSelfCorrectByOverlap();

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

		size_t SelectFreqsOfrange(const size_t LowerBound, const size_t UpperBound, leafList& newLeaves);

		std::pair<size_t,size_t> alnscore;
    private:

		//
		// Functions
		//
			void initialRootNode(const std::string& beginningkmer);
			void buildOverlapbyFMindex(IntervalTree<size_t>& fwdIntervalTree,IntervalTree<size_t>& rvcIntervalTree,const int& overlapSize);

			void extendLeaves(leafList& newLeaves);
			void attempToExtend(leafList& newLeaves, const bool isSuccessToReduce);
			void updateLeaves(leafList& newLeaves,extArray& extensions,leafInfo& leaf,size_t currLeavesNum);

			void refineSAInterval(leafList& leaves, const size_t newKmerSize);

			int findTheBestPath(const SAIntervalNodeResultVector& results, FMWalkResult2& FMWResult);

			extArray getFMIndexExtensions(const leafInfo& currLeaf,const bool printDebugInfo);

			// prone the leaves without seeds in proximity
				bool PrunedBySeedSupport(leafList& newLeaves);
			//Check if need reduce kmer size
				bool isInsufficientFreqs(leafList& newLeaves);

			// Check if the leaves reach $
				bool isTerminated(SAIntervalNodeResultVector& results);

			bool isOverlapAcceptable(SAIOverlapNode3* currNode);
			bool isSupportedByNewSeed(SAIOverlapNode3* currNode, size_t smallSeedIdx, size_t largeSeedIdx);
			match_t ismatchedbykmer(BWTInterval currFwdInterval,BWTInterval currRvcInterval);
			double computeErrorRate(SAIOverlapNode3* currNode);

			void printDebugData(debugPerExtArray& currDebugData);
			void printErrorRate(leafList& currLeaves);

		//
		// Data
		//
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
			size_t m_step_number;

		size_t m_maxIndelSize;
		double* freqsOfKmerSize;

		// Optional parameters
			size_t m_maxfreqs;

		std::string m_query;
		size_t m_maxLength;
		size_t m_minLength;
		std::vector<BWTInterval> m_fwdTerminatedInterval;   //in rBWT
		std::vector<BWTInterval> m_rvcTerminatedInterval;   //in BWT

		leafList m_leaves;
		SAIOverlapNode3* m_pRootNode;
		SONode3PtrList m_RootNodes;

		size_t m_currentLength;
		size_t m_currentKmerSize;

		IntervalTree<size_t> m_fwdIntervalTree;
		IntervalTree<size_t> m_rvcIntervalTree;

		IntervalTree<size_t> m_fwdIntervalTree2;
		IntervalTree<size_t> m_rvcIntervalTree2;

		size_t RemainedMaxLength;

};

#endif