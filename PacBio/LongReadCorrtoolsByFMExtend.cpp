#include "BWT.h"
#include "BWTAlgorithms.h"
#include "SAINode.h"
#include "HashtableSearch.h"
#include "LongReadCorrtoolsByFMExtend.h"

// simpleSeedFeature
	simpleSeedFeature::simpleSeedFeature
						(
							const SeedFeature& seed,
							bool isHead
						):
							seq     (seed.seedStr),
							end     (seed.seedEndPos),
							isRepeat(seed.isRepeat)
						{
										if (isHead)
											start = seed.seedEndPos - seed.seedLen+1;
										else
											start = seed.seedStartPos;
						}

// seedPair
	seedPair::seedPair  (
							const SeedFeature& source,
							const SeedFeature& target,
							bool isPosStrand
						):
							source(source,true),
							target(target,false),
							isPosStrand(isPosStrand) {};

	void seedPair::reverseSeed()
						{
							std::swap(source,target);
							source.seq = reverseComplement(source.seq);
							target.seq = reverseComplement(target.seq);
							isPosStrand = !isPosStrand;
						}

	void seedPair::reduceSourceBy(size_t length)
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

// FMidx_t
	FMidx_t::FMidx_t(
						const std::string& s,
						const BWTInterval& fwdInterval,
						const BWTInterval& rvcInterval
					):
						SearchLetters(s),
						fwdInterval(fwdInterval),
						rvcInterval(rvcInterval),
						kmerFrequency(fwdInterval.size() + rvcInterval.size())
					{};

	FMidx_t::FMidx_t(
						const char c,
						const BWTInterval& fwdInterval,
						const BWTInterval& rvcInterval
					):
						SearchLetters(std::string(1,c)),
						fwdInterval(fwdInterval),
						rvcInterval(rvcInterval),
						kmerFrequency(fwdInterval.size() + rvcInterval.size())
					{};

	void FMidx_t::setInterval(const BWTInterval& fwdInterval,const BWTInterval& rvcInterval)
			{
				this -> fwdInterval   = fwdInterval;
				this -> rvcInterval   = rvcInterval;
				this -> kmerFrequency = fwdInterval.size() + rvcInterval.size();
			}

	BWTInterval FMidx_t::getFwdInterval()
		{ return fwdInterval; }

	BWTInterval FMidx_t::getRvcInterval()
		{ return rvcInterval; }

	kmerFreq_t FMidx_t::getKmerFrequency()
		{ return kmerFrequency; }

	char FMidx_t::getLetter()
		{ return SearchLetters.front(); }

	bool FMidx_t::isVaildIntervals()
		{ return (fwdInterval.isValid() || rvcInterval.isValid()); }
// pathInfo_t
	pathInfo_t::pathInfo_t
						(const std::string& beginStr):
							localErrorRate (0),
							globalErrorRate(0)
						{ setFirstBaseInfo(beginStr); };

	pathInfo_t::pathInfo_t
						(
							const pathInfo_t& lastNode,
							const FMidx_t& currExt
						):
							localErrorRate (0),
							globalErrorRate(0)
						{  setNextBaseInfo(lastNode, currExt); }

	void pathInfo_t::setFirstBaseInfo
						(
							const std::string& beginStr
						)
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
								for(auto    reverseIdx = beginStr.crbegin();
									reverseIdx != beginStr.crend(); ++reverseIdx)
								{
									char suffixLetter = (*reverseIdx);
									if (reverseIdx == beginStr.crbegin())
										lastBase =  suffixLetter;
									if( lastBase == suffixLetter )
										Homopolymer++;
									else
										break;
								}
						}

	void pathInfo_t::setNextBaseInfo
						(
							const pathInfo_t& lastNode,
							const FMidx_t& currExt
						)
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
						}

	void pathInfo_t::setErrorRate
						(
							const double localErrorRate,
							const double globalErrorRate
						)
						{
							this -> localErrorRate  = localErrorRate;
							this -> globalErrorRate = globalErrorRate;
						}

// leafElement_t
	leafElement_t::leafElement_t
						(
							SAIOverlapNode3* node,
							const size_t     start
						):
							leaf    (node),
							birthday(start) {};

// SAIOverlapNode3
	SAIOverlapNode3::SAIOverlapNode3
						(
							const std::string* pQuery,
							SAIOverlapNode3* parent
						):
							SAINode(pQuery,parent),
							lastLeafID     (0)  ,
							lastOverlapLen (0)  , currOverlapLen   (0), queryOverlapLen(0),
							lastSeedIdx    (0)  , lastSeedIdxOffset(0),
							numRedeemSeed  (0)  , totalSeeds       (0), numOfErrors    (0)
						{};

	SAIOverlapNode3::~SAIOverlapNode3()
						{
							// Delete children
								for(auto iter = m_children3.begin(); iter != m_children3.end(); ++iter)
									delete *iter;
						};
	// Add a child node to this node with the given label and return a pointer to the created node
	SAIOverlapNode3* SAIOverlapNode3::createChild(const std::string& label)
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

	// get the ancestorNode corresponding the leaf at a certain location.
	SAIOverlapNode3* SAIOverlapNode3::findAncestor(const size_t location)
						{
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

								return ancestorNode;
						}
	// Set root
	void SAIOverlapNode3::setRoot
						(
							const std::string& beginStr,
							const size_t& beginLen,
							const size_t& matchSize,
							const BWT* m_pBWT,
							const BWT* m_pRBWT
						)
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
	void SAIOverlapNode3::setNode
						(
							FMidx_t& extension,
							const size_t currLeavesNum,
							SAIOverlapNode3* refNode,
							const size_t m_step_number
						)
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

	pathInfo_t& SAIOverlapNode3::getParentInfo
						(
							pathPack_t& pathPack,
							const size_t location
						)
						{
							if (location == (size_t) -1)
								return pathPack.back().at(this);

							// get the ancestorNode corresponding the leaf at a certain location.
								SAIOverlapNode3* ancestorNode = findAncestor(location);

							// get the goal-parent information
								return pathPack.at(location).at(ancestorNode);
						}

// dominantBase
	dominantBase::dominantBase
						(const std::string& goalBase):
							identicalBase(goalBase.front()) {};

	void dominantBase::setFreq(FMidx_t currExt)
			{
				kmerFreq_t currFreq = currExt.getKmerFrequency();
				char       currBase = currExt.getLetter();

				if( currFreq >  freqs.at(0)
				|| (currFreq == freqs.at(0) && currBase == identicalBase))
				{
					std::swap(freqs.at(0), currFreq);
					std::swap(bases.at(0), currBase);
				}
				if( currFreq >  freqs.at(1)
				|| (currFreq == freqs.at(1) && currBase == identicalBase))
				{
					std::swap(freqs.at(1), currFreq);
					std::swap(bases.at(1), currBase);
				}
			}

	bool dominantBase::isDominant()
		{ return (bases.at(0) == identicalBase || bases.at(1) == identicalBase); }

	kmerFreq_t dominantBase::getMaxFreq()
		{ return freqs.front(); }

// forceExtInfo
	void forceExtInfo::setSpecPath
						(
							const bool& isSpecifiedPath,
							const std::string& specifiedPath
						)
						{
							this -> isSpecifiedPath = isSpecifiedPath;
							this -> specifiedPath   = specifiedPath;
						};

	bool forceExtInfo::isEmpty()
		{ return specifiedPath.empty(); }

	bool forceExtInfo::isMatchBase(const std::string& currBase, const size_t position)
		{ return (currBase.back() == specifiedPath.at(position));}

	bool forceExtInfo::isMatchBase(const char currBase, const size_t position)
		{ return (currBase == specifiedPath.at(position));}

	bool forceExtInfo::isIllegalLen(size_t position)
		{ return (position >= specifiedPath.length());}

	void forceExtInfo::RCSeq()
		{
			if (!isEmpty())
				specifiedPath = reverseComplement(specifiedPath);
		}

	std::string forceExtInfo::getSpecifiedPath()
		{ return specifiedPath; }

// debugExtInfo
	debugExtInfo::debugExtInfo
						(
							bool isDebug,
							std::ostream* debug_file    ,
							std::ostream* debug_finalSeq,
							std::string readID,
							forceExtInfo ref
						):
							isDebug(isDebug),
							debug_file(debug_file),
							debug_finalSeq(debug_finalSeq),
							readID(readID),
							ref(ref) {};
// debugPerExtInfo
	debugPerExtInfo::debugPerExtInfo
						(
							FMidx_t currExt,
							bool    isMatched
						):
							isMatched(isMatched)
						{
							(this  -> kmerFreq ) = currExt.getKmerFrequency();
							(this -> errorType ).reserve(4);
							(this  -> kmerRatio) = 0;
						}

		void debugPerExtInfo::setErrorData
						(
							const double       kmerRatio,
							const std::string& errorString
						)
						{
							(this -> kmerRatio) = kmerRatio;
							(this -> errorType) = "X" + errorString;
						}
		void debugPerExtInfo::extSuccessSet()
						{
							(this -> errorType).front() = 'O';
						}