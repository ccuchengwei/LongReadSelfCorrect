#ifndef LONGREADCORRTOOLS_H
#define LONGREADCORRTOOLS_H

#include "SeedFeature.h"
#include <unordered_map>
#include <list>

// Define int type.
	typedef int64_t kmerFreq_t;

// Define class
	class FMWalkResult2;
	class FMextendParameters;

	class simpleSeedFeature;
	class seedPair;

	class FMidx_t;
	class pathInfo_t;
	class leafElement_t;
	class SAIOverlapNode3;

	class dominantBase;


	class forceExtInfo;
	class debugExtInfo;
	class debugPerExtInfo;

// Define class group.
	typedef std::vector<FMidx_t> extArray;

	typedef std::unordered_map<SAIOverlapNode3*, pathInfo_t> leafTable_t;
	typedef std::vector<leafTable_t>     pathPack_t;
	typedef std::vector<leafElement_t>   leafAncestry_t;
	typedef std::list<SAIOverlapNode3*>  SONode3PtrList;

	typedef std::vector<debugPerExtInfo> debugPerExtArray;

class FMWalkResult2
{
	public:
		std::string mergedSeq;
		int alnScore;
		double kmerFreq;
};

class FMextendParameters
{
	public:
		// Functions
			FMextendParameters  (
									BWTIndexSet  indices,
									int      idmerLength,
									int        maxLeaves,
									int    minKmerLength,
									size_t    PBcoverage,
									double ErrorRate
								):
									indices(indices),
									idmerLength(idmerLength),
									maxLeaves(maxLeaves),
									minKmerLength(minKmerLength),
									PBcoverage(PBcoverage),
									ErrorRate(ErrorRate)
								{};

			FMextendParameters(){};

		// Data
			BWTIndexSet indices;
			int idmerLength;
			int maxLeaves;
			int minKmerLength;
			size_t PBcoverage;
			double ErrorRate;

};

class simpleSeedFeature
{
	public:
		// Functions
			simpleSeedFeature(const SeedFeature& seed, bool isHead);

		// Data
			std::string seq;
			size_t start;
			size_t end;
			bool   isRepeat;
};

class seedPair
{
	public:
		// Functions
			seedPair(
						const SeedFeature& source,
						const SeedFeature& target,
						bool isPosStrand = true
					);

			void reverseSeed();
			void reduceSourceBy(size_t length);

		// Data
			simpleSeedFeature  source;
			simpleSeedFeature  target;
			bool    isPosStrand;
};

class FMidx_t
{
	public:
		// Functions
			FMidx_t (
						const std::string& s,
						const BWTInterval& fwdInterval,
						const BWTInterval& rvcInterval
					);

			FMidx_t (
						const char c,
						const BWTInterval& fwdInterval,
						const BWTInterval& rvcInterval
					);

			void setInterval(const BWTInterval& fwdInterval,const BWTInterval& rvcInterval);
			BWTInterval getFwdInterval();
			BWTInterval getRvcInterval();
			kmerFreq_t getKmerFrequency();
			char getLetter();
			bool isVaildIntervals();

		// Data
			const std::string SearchLetters;

	private:
		BWTInterval fwdInterval;
		BWTInterval rvcInterval;
		kmerFreq_t kmerFrequency;

};

class pathInfo_t
{
	public:
		// Functions
			pathInfo_t  () = default;
			pathInfo_t  (const std::string& beginStr);
			pathInfo_t  (const pathInfo_t& lastNode, const FMidx_t& currExt);
			void setFirstBaseInfo (const std::string& beginStr);
			void setNextBaseInfo  (const pathInfo_t& lastNode, const FMidx_t& currExt);
			void setErrorRate(const double localErrorRate, const double globalErrorRate);

		// Data
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
		// Functions
			leafElement_t   (
								SAIOverlapNode3* node,
								const size_t     start
							);
		// Data
			SAIOverlapNode3* leaf;
			size_t           birthday;
};

class SAIOverlapNode3 : public SAINode
{
	public:
		// Functions
				SAIOverlapNode3(const std::string* pQuery, SAIOverlapNode3* parent);
				~SAIOverlapNode3();
				SAIOverlapNode3* createChild (const std::string& label);
				SAIOverlapNode3* findAncestor(const size_t location);

				void setRoot(
								const std::string& beginStr,
								const size_t& beginLen,
								const size_t& matchSize,
								const BWT* m_pBWT,
								const BWT* m_pRBWT
							);

				void setNode(
								FMidx_t& extension,
								const size_t currLeavesNum,
								SAIOverlapNode3* refNode,
								const size_t m_step_number = 0
							);

				pathInfo_t& getParentInfo
							(
								pathPack_t& pathPack,
								const size_t location = (size_t) -1
							);
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

class dominantBase
{
	public:
		// Functions
			dominantBase(const std::string& goalBase);
			void setFreq(FMidx_t currExt);
			bool isDominant();
			kmerFreq_t getMaxFreq();

		// Data
			std::vector<kmerFreq_t>  freqs = {0,0};
			std::string              bases = {0,0};
			char identicalBase = 0;
};

class forceExtInfo
{
	public:
		// Functions
			forceExtInfo():
							isSpecifiedPath(false),
							specifiedPath(""){};

			void setSpecPath(
								const bool& isSpecifiedPath,
								const std::string& specifiedPath
							);

			bool isEmpty();
			bool isMatchBase(const std::string& currBase, const size_t position);
			bool isMatchBase(const char currBase        , const size_t position);
			bool isIllegalLen(size_t position);
			void RCSeq();
		std::string getSpecifiedPath();

		// Data
			bool isSpecifiedPath;
	private:
		std::string specifiedPath;

};

class debugExtInfo
{
	public:
		debugExtInfo(
						bool isDebug = false,
						std::ostream* debug_file     = nullptr,
						std::ostream* debug_finalSeq = nullptr,
						std::string readID = "",
						forceExtInfo ref = forceExtInfo()
					);

		const bool isDebug;
		std::ostream* debug_file;
		std::ostream* debug_finalSeq;
		std::string   readID;
		forceExtInfo  ref;

};

class debugPerExtInfo
{
	public:
		// Functions
			debugPerExtInfo (
								FMidx_t currExt,
								bool    isMatched = 0
							);

			void setErrorData   (
									const double       kmerRatio,
									const std::string& errorString
								);

			void extSuccessSet();

		// Data
			std::string errorType;
			kmerFreq_t  kmerFreq;
			double      kmerRatio;
			bool        isMatched;
};

#endif
