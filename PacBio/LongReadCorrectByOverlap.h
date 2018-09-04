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

#include <unordered_map>
#include <list>
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SAINode.h"
#include "IntervalTree.h"
#include "HashtableSearch.h"
#include "SeedFeature.h"


// Define int type.
	typedef int64_t kmerFreq_t;

// Define class
	class pathInfo_t;
	class leafElement_t;
	class SAIOverlapNode3;

// Define class group.
	typedef std::unordered_map<SAIOverlapNode3*, pathInfo_t> leafTable_t;
	typedef std::vector<leafTable_t>    pathPack_t;
	typedef std::vector<leafElement_t>  leafAncestry_t;
	typedef std::list<SAIOverlapNode3*> SONode3PtrList;

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

class pathInfo_t
{
	public:
		// function
			pathInfo_t  () = default;

			pathInfo_t  (const std::string& beginStr):
								localErrorRate (0),
								globalErrorRate(0)
						{

							// Set GC Number
								GCNumber = 0;
								for(auto currBase : beginStr)
								{
									if (currBase == 'C' || currBase == 'G')
										GCNumber++;
								}
							// Set tailLetter and its counter
								Homopolymer = 0;
								for(auto reverseIdx = beginStr.crbegin(); reverseIdx != beginStr.crend(); ++reverseIdx)
								{
									char suffixLetter = (*reverseIdx);
									if (reverseIdx == beginStr.crbegin())
										lastBase =  suffixLetter;
									if( lastBase == suffixLetter )
										Homopolymer++;
									else
										break;
								}
						};

				pathInfo_t  (const pathInfo_t& lastNode, const FMidx_t& currExt):
								localErrorRate (0),
								globalErrorRate(0)
						{
							char currBase = currExt.SearchLetters.back();

							// Set GC Number
								GCNumber = lastNode.GCNumber;
								if (currBase == 'C' || currBase == 'G')
									GCNumber++;

							// Set tailLetter and its counter
								lastBase = currBase;
								if (lastNode.lastBase == currBase)
									Homopolymer = (lastNode.Homopolymer) + 1;
								else
									Homopolymer = 1;
						};

				void setErrorRate(const double localErrorRate, const double globalErrorRate)
					{
						this -> localErrorRate  = localErrorRate;
						this -> globalErrorRate = globalErrorRate;
					}
		// data
			// last Base
				char lastBase;
			// Error Rates
				double localErrorRate;
				double globalErrorRate;
			// the parameters to detect the regions where tends to errors
				size_t GCNumber;
				size_t Homopolymer;

};

class leafElement_t
{
	public:
		// function
			leafElement_t   (
								SAIOverlapNode3* node,
								const size_t     start
							):
								leaf    (node),
								birthday(start)
								{};
		// data
			SAIOverlapNode3* leaf;
			size_t           birthday;
};

class SAIOverlapNode3 : public SAINode
{
	public:
		// Functions
				SAIOverlapNode3(const std::string* pQuery, SAIOverlapNode3* parent):
					SAINode(pQuery,parent),
					lastLeafID     (0)  ,
					lastOverlapLen (0)  , currOverlapLen   (0), queryOverlapLen(0),
					lastSeedIdx    (0)  , lastSeedIdxOffset(0),
					numRedeemSeed  (0)  , totalSeeds       (0), numOfErrors    (0)
					{};

				~SAIOverlapNode3()
				{
					// Delete children
					for(SONode3PtrList::iterator iter = m_children3.begin(); iter != m_children3.end(); ++iter)
						delete *iter;
				};

			// Add a child node to this node with the given label
			// Returns a pointer to the created node
				SAIOverlapNode3* createChild(const std::string& label)
				{
					SAIOverlapNode3* pAdded = new SAIOverlapNode3(m_pQuery, this);
					pAdded->extend(label);

					// still lack of a copy constructor
						pAdded->lastOverlapLen        = lastOverlapLen;
						pAdded->currOverlapLen        = currOverlapLen;
						pAdded->queryOverlapLen       = queryOverlapLen;
						pAdded->initSeedIdx           = initSeedIdx;
						pAdded->lastSeedIdx           = lastSeedIdx;
						pAdded->lastSeedIdxOffset     = lastSeedIdxOffset;
						pAdded->numRedeemSeed         = numRedeemSeed;
						pAdded->totalSeeds            = totalSeeds;
						pAdded->numOfErrors           = numOfErrors;
						pAdded->Ancestor              = Ancestor;
						pAdded->resultindex           = resultindex;

						// pAdded->currkmersize = this->currkmersize;
						m_children3.push_back(pAdded);

					return pAdded;
				}

			// Set root
				void setRoot(const std::string& beginStr, const size_t& beginLen, const size_t& matchSize, const BWT* m_pBWT, const BWT* m_pRBWT)
				{
					lastLeafID  = 1;
					// Copy the intervals
						computeInitial(beginStr);
						fwdInterval = BWTAlgorithms::findInterval(m_pRBWT, reverse(beginStr));
						rvcInterval = BWTAlgorithms::findInterval(m_pBWT , reverseComplement(beginStr));
						addKmerCount( fwdInterval.size() + rvcInterval.size() );

					queryOverlapLen =
							currOverlapLen =
							lastOverlapLen = beginLen;

					lastSeedIdx = initSeedIdx = beginLen - matchSize;

					numRedeemSeed = 0;
					totalSeeds = beginLen - matchSize + 1;

					Ancestor.emplace_back(this, 0);

				}

			// Set children
				void setNode(FMidx_t& extension, const size_t currLeavesNum, SAIOverlapNode3* refNode, const size_t m_step_number = 0)
				{
					// Set currNode
						// Set lastLeafID
							lastLeafID = currLeavesNum;
						// Copy the intervals
							fwdInterval = extension.getFwdInterval();
							rvcInterval = extension.getRvcInterval();
							addKmerCount( extension.getKmerFrequency() );
						// currOverlapLen/queryOverlapLen always increase wrt each extension
						// in order to know the approximate real-time matched length for terminal/containment processing
							currOverlapLen++;
							queryOverlapLen++;
						// Set the Ancestors
							if (this != refNode)
								Ancestor.emplace_back(this, m_step_number);
				}

				pathInfo_t& getParentInfo(pathPack_t& pathPack, const size_t location = (size_t) -1)
				{
					if (location == (size_t) -1)
						return pathPack.back().at(this);

					// get the ancestorNode corresponding the leaf at a certain location.
						SAIOverlapNode3* ancestorNode   = Ancestor.front().leaf;

						// binary search to find ancestorNode.
						size_t start  = 0;
						size_t end    = Ancestor.size() -1;
						size_t middle = (start + end)/2;
						while(start <= end)
						{
							middle = (start + end)/2;

							if (Ancestor.at(middle).birthday <= location)
							{
								ancestorNode = Ancestor.at(middle).leaf;
								start = middle + 1;
							}
							else
								end = middle - 1;
						}

					// get the goal-parent information
						return pathPack.at(location).at(ancestorNode);
				}
		// Variables
			// Leaf ID
				size_t lastLeafID;
			// BWTIntervals
				BWTInterval fwdInterval;
				BWTInterval rvcInterval;
			// Match Variables
				// lengths
					// last overlap length when matching last seed
						size_t lastOverlapLen;
					// current overlap length on the subject increases wrt each FM-index extension
						size_t currOverlapLen;
					// current overlap length on the query
						size_t queryOverlapLen;
				// locations
					// index of the init seed
						int initSeedIdx;
					// last matched seed index
						size_t lastSeedIdx;
					// error seed begin idx
						// size_t errorSeedBeginIdx;
					// index offset to the center
						int lastSeedIdxOffset;
				// match
					// number of redeem seeds
						double numRedeemSeed;
					// number of real matches
						size_t totalSeeds;
					// number of SNPs or indels
						size_t numOfErrors;
				// the total path in the current leaf.
					leafAncestry_t Ancestor;

			// Result
				// index of the result and index of the matchpoint
				std::pair <int,int> resultindex = std::make_pair(-1,-1);
			// Others
				// size_t currkmersize;
			// Children
				SONode3PtrList m_children3;
};

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

struct debugPerExtInfo
{
	debugPerExtInfo (
						FMidx_t currExt,
						bool    isMatched = 0
					):
						isMatched(isMatched)
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
	bool        isMatched;
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