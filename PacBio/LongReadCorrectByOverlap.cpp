///----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------
//
// LongReadSelfCorrectByOverlap - A fast overlap filter using locality-sensitive backward search.
//
//
#include <iomanip>
#include "LongReadCorrectByOverlap.h"
#include "KmerFeature.h"
#include "BWTAlgorithms.h"
#include "stdaln.h"
#include "LongReadOverlap.h"

// Class: SAIOverlapTree
LongReadSelfCorrectByOverlap::LongReadSelfCorrectByOverlap
				(
					const seedPair& extSeeds,
					const std::string& strBetweenSrcTarget,
					int disBetweenSrcTarget,
					size_t initkmersize,
					size_t maxOverlap,
					const FMextendParameters params,
					size_t min_SA_threshold,
					const debugExtInfo debug,
					double errorRate,
					size_t repeatFreq,
					size_t localSimilarlykmerSize
				):
					m_extSeeds(extSeeds),
					m_sourceSeed(extSeeds.source.seq),
					m_strBetweenSrcTarget(strBetweenSrcTarget),
					m_targetSeed(extSeeds.target.seq),
					m_disBetweenSrcTarget(disBetweenSrcTarget),
					m_initkmersize(initkmersize),
					m_minOverlap(params.minKmerLength),
					m_maxOverlap(maxOverlap),
					m_pBWT(params.indices.pBWT),
					m_pRBWT(params.indices.pRBWT),
					m_PBcoverage(params.PBcoverage),
					m_min_SA_threshold(min_SA_threshold),
					m_errorRate(errorRate),
					m_maxLeaves(params.maxLeaves),
					m_seedSize(params.idmerLength),
					m_repeatFreq(repeatFreq),
					m_localSimilarlykmerSize(localSimilarlykmerSize),
					m_PacBioErrorRate(params.ErrorRate),
					m_Debug(debug)
{
	m_extSeeds.reduceSourceBy(m_sourceSeed.length()-m_initkmersize);

	//if distance < 100 ,use const indel size
		if (m_disBetweenSrcTarget > 100)
			m_maxIndelSize  =  m_disBetweenSrcTarget * 0.2;
		else
			m_maxIndelSize  =  20;

	//initialRootNode
		initialRootNode(m_extSeeds.source.seq);

	// push new node into roots and leaves vector
		m_RootNodes.push_back(m_pRootNode);
		m_leaves.emplace_back(m_pRootNode,1);

	//frequencies of correspond k
		freqsOfKmerSize = new double[100 + 1]{0};
		for(int i = m_minOverlap ; i <= 100 ; i++)
			freqsOfKmerSize[i] = pow(1 - m_PacBioErrorRate, i) * m_PBcoverage;
//*
	if(m_Debug.isDebug)
	{
		std::cout   << m_Debug.caseNum << "\tBE: "
					<< m_extSeeds.source.seq << " "
					<< m_extSeeds.target.seq << " "
					<< m_pRootNode->fwdInterval.size() + m_pRootNode->rvcInterval.size() << "|"
					<< disBetweenSrcTarget <<"\n";
	}
//*/
	// PacBio reads are longer than real length due to insertions
		m_maxLength =          (1.2*(m_disBetweenSrcTarget+10))+2*m_initkmersize     ;
		m_minLength = std::max((0.8*(m_disBetweenSrcTarget-20))+2*m_initkmersize,0.0);

	// initialize the ending SA intervals with kmer length = m_minOverlap
		for(size_t i =0 ;i <= m_targetSeed.length()-m_minOverlap; i++)
		// for(size_t i =0 ;i <= 0; i++)
		{
			std::string endingkmer = m_targetSeed.substr(i, m_minOverlap);

			m_fwdTerminatedInterval.push_back(BWTAlgorithms::findInterval(m_pRBWT, reverse(endingkmer)));
			m_rvcTerminatedInterval.push_back(BWTAlgorithms::findInterval(m_pBWT, reverseComplement(endingkmer)));
		}
    // build overlap tree
		m_query = m_extSeeds.source.seq + m_strBetweenSrcTarget + m_targetSeed;
		// build overlap tree to determine the error rate
			buildOverlapbyFMindex(m_fwdIntervalTree ,m_rvcIntervalTree ,m_seedSize);
		// build overlap tree to match 5-mer
			buildOverlapbyFMindex(m_fwdIntervalTree2,m_rvcIntervalTree2,5);

}

LongReadSelfCorrectByOverlap::~LongReadSelfCorrectByOverlap()
{
	for (auto iter = m_RootNodes.begin(); iter != m_RootNodes.end(); ++iter)
		delete *iter;

	m_RootNodes.clear();

	delete[] freqsOfKmerSize;
}

// Initialize the root node
void LongReadSelfCorrectByOverlap::initialRootNode(const std::string& beginningkmer)
{
	// create one root node
		m_pRootNode = new SAIOverlapNode3(&m_sourceSeed, NULL);

	// store initial str of root
		m_pRootNode->computeInitial(beginningkmer);
		m_pRootNode->fwdInterval = BWTAlgorithms::findInterval(m_pRBWT, reverse(beginningkmer));
		m_pRootNode->rvcInterval = BWTAlgorithms::findInterval(m_pBWT , reverseComplement(beginningkmer));
		m_pRootNode->lastOverlapLen = m_currentLength = m_pRootNode->currOverlapLen = m_pRootNode->queryOverlapLen = m_currentKmerSize = m_initkmersize;
		m_pRootNode->lastSeedIdx = m_pRootNode->initSeedIdx = m_initkmersize - m_seedSize;
		m_pRootNode->totalSeeds = m_initkmersize - m_seedSize + 1;
		m_pRootNode->numRedeemSeed = 0;
		m_pRootNode->LocalErrorRateRecord.push_back(0);
		m_pRootNode->GlobalErrorRateRecord.push_back(0);
		m_maxfreqs = m_pRootNode->fwdInterval.size() + m_pRootNode->rvcInterval.size();
}

//Build the overlap tree
void LongReadSelfCorrectByOverlap::buildOverlapbyFMindex(IntervalTree<size_t>& fwdIntervalTree,IntervalTree<size_t>& rvcIntervalTree,const int& overlapSize)
{
	// put SA intervals into fwdIntervals and rvcIntervals cache
		std::vector< TreeInterval<size_t> > fwdIntervals;
			fwdIntervals.reserve(m_query.length()-overlapSize+1);

		std::vector< TreeInterval<size_t> > rvcIntervals;
			rvcIntervals.reserve(m_query.length()-overlapSize+1);

	// Build the Overlap Tree
		for(int i = 0; i <= (int)m_query.length()-(int)overlapSize ; i++)//build overlap tree
		{

			std::string seedStr = m_query.substr(i, overlapSize);

			BWTInterval bi;
			bi = BWTAlgorithms::findInterval( m_pRBWT, reverse(seedStr) );
			if(bi.isValid())
				fwdIntervals.emplace_back( bi.lower, bi.upper, i );
			bi = BWTAlgorithms::findInterval( m_pBWT, reverseComplement(seedStr) );
			if(bi.isValid())
				rvcIntervals.emplace_back( bi.lower, bi.upper, i );
		}
		fwdIntervalTree = IntervalTree<size_t>(fwdIntervals);
		rvcIntervalTree = IntervalTree<size_t>(rvcIntervals);
}

//On success return the length of merged string
int LongReadSelfCorrectByOverlap::extendOverlap(FMWalkResult2& FMWResult)
{
	SAIntervalNodeResultVector results;
	m_step_number = 1;

	//Overlap extension via FM-index walk
	while(!m_leaves.empty() && m_leaves.size() <= m_maxLeaves && m_currentLength <= m_maxLength)
	{
/*
		if(m_Debug.isDebug)
		{
			std::cout   << m_Debug.caseNum << "\t" << m_step_number
						<< " Leaves number for extension:" << m_leaves.size() << std::endl;
			std::cout << "\t(Rel) Leaves number to trim branch:" << m_leaves.size() << std::endl;
			std::cout << "\t";
			printErrorRate(m_leaves);
		}
/*/
		// ACGT-extend the leaf nodes via updating existing SA interval
			leafList newLeaves;
			extendLeaves(newLeaves);
/*

		if(m_Debug.isDebug)
		{
			std::cout << "\t(Rel) Leaves number after trimming branch: " << m_leaves.size() << std::endl;
			std::cout << "\t(Abs) Leaves number to trim branch:" << newLeaves.size() << std::endl;
			std::cout << "\t";
			printErrorRate(newLeaves);
		}
//*/
		//use overlap tree to trim branch
			PrunedBySeedSupport(newLeaves);
/*

		if(m_Debug.isDebug)
		{
			std::cout << "\t(Abs) Leaves number after trimming branch: " << newLeaves.size() << std::endl;
			std::cout << "\tCurrent Length:" << m_currentLength << " (" << m_minLength << "," << m_maxLength << ")" << std::endl;
		}
//*/
		//update leaves
			m_leaves.clear();
			m_leaves = newLeaves;
/*
		if(m_Debug.isDebug)
			std::cout << "----" << std::endl;
//*/
		if(m_currentLength >= m_minLength)
			isTerminated(results);

		m_step_number++;

		if (m_Debug.ref.isSpecifiedPath)
		{
			if (m_Debug.ref.isIllegalLen(m_step_number-1))
			{
				m_leaves.clear();
				break;
			}
		}

	}

	// reach the terminal kmer
	if(results.size() > 0)
		return findTheBestPath(results, FMWResult); // find the path with maximum match percent or kmer coverage
/*
	if(m_Debug.isDebug)
			std::cout << "\tERROR\t" << std::endl;
//*/

	// Did not reach the terminal kmer
	if(m_leaves.empty())	//high error
		return -1;
	else if(m_currentLength > m_maxLength)	//exceed search depth
		return -2;
	else if(m_leaves.size() > m_maxLeaves)	//too much repeats
		return -3;
	else
		return -4;
}

// find the path where is the min error rate.
int LongReadSelfCorrectByOverlap::findTheBestPath(const SAIntervalNodeResultVector& results, FMWalkResult2& FMWResult)
{
	double minErrorRate = 1;

	for (size_t i = 0 ; i < results.size() ;i++)
	{
		const std::string& candidateSeq = results[i].thread;
//*
		if(m_Debug.isDebug)
			{
				std::cout << "Final Seqs:" << results[i].thread << "\tError Rate:" << results[i].errorRate << std::endl;
			}
//*/
		if(results[i].errorRate < minErrorRate )
		{
			minErrorRate = results[i].errorRate;
			FMWResult.mergedSeq = candidateSeq;
			minTotalcount = results[i].SAIntervalSize;
		}
	}

	if(FMWResult.mergedSeq.length() != 0)
		return 1;
	return -4;
}

// Extend all Leaves by FM-index
void LongReadSelfCorrectByOverlap::extendLeaves(leafList& newLeaves)
{
	//resize if length too long
	if(m_currentKmerSize > m_maxOverlap)
		refineSAInterval(m_leaves, m_maxOverlap);

	attempToExtend(newLeaves,true);

	if(newLeaves.empty() ) //level 1 reduce size
	{
		size_t LowerBound = std::max(m_currentKmerSize - 2, m_minOverlap);
		size_t ReduceSize = SelectFreqsOfrange(LowerBound,m_currentKmerSize,m_leaves);
		bool isSuccessToReduce = m_currentKmerSize != ReduceSize;
		refineSAInterval(m_leaves, ReduceSize);

		attempToExtend(newLeaves,isSuccessToReduce);

		if( newLeaves.empty() )//level 2 reduce threshold
		{
			m_min_SA_threshold--;
			attempToExtend(newLeaves,false);
			m_min_SA_threshold++;
		}
	}

	//extension succeed
	if(!newLeaves.empty())
	{
		m_currentLength++;
		m_currentKmerSize++;
		if( isInsufficientFreqs(newLeaves) )// if frequency are low , relax it
		{
			size_t LowerBound = std::max(m_currentKmerSize - 2, m_minOverlap);
			size_t ReduceSize = SelectFreqsOfrange(LowerBound,m_currentKmerSize,newLeaves);
			refineSAInterval(newLeaves,ReduceSize);
		}

	}

}

// Determine the size of the reduced k-mer in leaves when the max freq of them approach to the expected one.
size_t LongReadSelfCorrectByOverlap::SelectFreqsOfrange(const size_t LowerBound, const size_t UpperBound, leafList& newLeaves)
{
	extArray maxKmerArray;
	int tempmaxfmfreqs = 0;

	for(auto& iter : newLeaves)
	{
		SAIOverlapNode3* leaf    = iter.leafNodePtr;
		std::string    maxKmer = leaf -> getSuffix(UpperBound);

		std::string startkmer  = maxKmer.substr(UpperBound - LowerBound); //  string of lower bound kmer size

		BWTInterval Fwdinterval = BWTAlgorithms::findInterval(m_pBWT, startkmer);
		BWTInterval Rvcinterval = BWTAlgorithms::findInterval(m_pRBWT, reverseComplement(reverse(startkmer)));

		maxKmerArray.emplace_back(maxKmer,Fwdinterval,Rvcinterval);
		FMidx& currKmer = maxKmerArray.back();

		if(currKmer.getKmerFrequency() > tempmaxfmfreqs ) //check interval size
			tempmaxfmfreqs = currKmer.getKmerFrequency();

	}

	if( tempmaxfmfreqs - (int)freqsOfKmerSize[LowerBound] < 5 ) return LowerBound;

	for(size_t i=1 ; i <= UpperBound - LowerBound; i++ )
	{
		tempmaxfmfreqs = 0;
		for(size_t j = 0; j < maxKmerArray.size(); j++)
		{
			std::string startkmer   = maxKmerArray.at(j).SearchLetters.substr(UpperBound - LowerBound - i);
			BWTInterval Fwdinterval = maxKmerArray.at(j).getFwdInterval();
			BWTInterval Rvcinterval = maxKmerArray.at(j).getRvcInterval();

			char b   = startkmer[0];//b:base
			char rcb = complement(b);
			BWTAlgorithms::updateInterval(Fwdinterval,  b,m_pBWT );
			BWTAlgorithms::updateInterval(Rvcinterval,rcb,m_pRBWT);

			maxKmerArray.at(j).setInterval(Fwdinterval,Rvcinterval);

			if(maxKmerArray.at(j).getKmerFrequency() > tempmaxfmfreqs) //check interval size
				tempmaxfmfreqs = maxKmerArray.at(j).getKmerFrequency();

		}

		if( tempmaxfmfreqs - (int)freqsOfKmerSize[LowerBound + i] < 5 ) return LowerBound + i ;
	}

	return UpperBound;
}

// Determine if frequencies in leaves are too low to perform the next extension.
bool LongReadSelfCorrectByOverlap::isInsufficientFreqs(leafList& newLeaves)
{
	size_t highfreqscount = 0;
	for(auto& iter : newLeaves)
	{
		int highfreqThreshold = m_PBcoverage > 60 ? (size_t)(m_PBcoverage/60)*3 : 3;

		if( iter.kmerFrequency > highfreqThreshold)
			highfreqscount++;
	}

	if(highfreqscount == 0)
		return true;
	else if (highfreqscount <= 2 && newLeaves.size()>=5 )
		return true;
	else if (highfreqscount <= 1 && newLeaves.size()>=3 )
		return true;
	return false;
}

// Refine SA intervals of each leave with a new k-mer
void LongReadSelfCorrectByOverlap::refineSAInterval(leafList& leaves, const size_t newKmerSize)
{
	for(auto& iter : leaves)
	{
		SAIOverlapNode3* leaf = iter.leafNodePtr;

		// reset the SA intervals using newKmerSize
			std::string reducedKmer = leaf -> getSuffix(newKmerSize);

			leaf -> fwdInterval = BWTAlgorithms::findInterval(m_pRBWT,           reverse(reducedKmer));
			leaf -> rvcInterval = BWTAlgorithms::findInterval(m_pBWT , reverseComplement(reducedKmer));
	}

	m_currentKmerSize = newKmerSize;
}

// Keep the leaves whose highly-correct rates are relative to the others.
// And attempt to extend those leaves.
void LongReadSelfCorrectByOverlap::attempToExtend(leafList& newLeaves, const bool isSuccessToReduce)
{
	double minimumErrorRate = 1;
	m_maxfreqs = 0;

	std::vector<size_t> frequencies;

	// Compute the min error rate
	for(auto& iter : m_leaves)
	{
		SAIOverlapNode3* leaf = iter.leafNodePtr;
		if( leaf->LocalErrorRateRecord.back() < minimumErrorRate)
			minimumErrorRate = leaf->LocalErrorRateRecord.back();
	}

	// Compute the errorRateDiff to trim leaves whose error rates relative to the others is high.
	leafList::iterator iter = m_leaves.begin();
	while(iter != m_leaves.end())
	{
		if (m_Debug.ref.isSpecifiedPath)
			break;

		SAIOverlapNode3* leaf = (*iter).leafNodePtr;
		double errorRateDiff  = (leaf->LocalErrorRateRecord.back()) - minimumErrorRate;
		if((errorRateDiff > 0.05 && m_currentLength > m_localSimilarlykmerSize/2)
		|| (errorRateDiff > 0.1  && m_currentLength > 15))
		{
			iter = m_leaves.erase(iter);
			continue;
		}
		++iter;
	}

	minTotalcount = 10000000;
	size_t currLeavesNum = 1;

	iter = m_leaves.begin();
	while(iter != m_leaves.end())
	{
		extArray extensions;
		int count = 0;
		SAIOverlapNode3* leaf = (*iter).leafNodePtr;
		bool isReduced = isSuccessToReduce;
		while(count < 2)
		{
			if	( count == 1
			&& !(leaf->LocalErrorRateRecord.back() == minimumErrorRate && m_leaves.size() > 1))
				break;

			if (m_Debug.isDebug && isReduced)
			{
				int kmer_freq = leaf->fwdInterval.size() + leaf->rvcInterval.size();

				char strand;
				if(m_extSeeds.isPosStrand)
					strand = '+';
				else
					strand = '-';

				*(m_Debug.debug_file)   << ">" << m_Debug.readID
										<< "|" << m_extSeeds.source.start
										<< "|" << m_extSeeds.source.end
										<< "|" << m_extSeeds.target.start
										<< "|" << m_extSeeds.target.end
										<< "|" << strand
										<< "&" << m_Debug.caseNum     << "|" << m_step_number
										<< "&" << m_leaves.size()
										<< "|" << currLeavesNum       << "|" << ((*iter).lastLeafID)
										<< "&" << m_currentKmerSize   << "|" <<  kmer_freq
										<< "&" << std::fixed          << std::setprecision(5)
										<< 100*(leaf ->LocalErrorRateRecord.back() ) << "|"
										<< 100*(leaf ->GlobalErrorRateRecord.back()) << "&";
			}

			extensions = getFMIndexExtensions(*iter,isReduced);

			if (m_Debug.isDebug && isReduced)
			{
				*(m_Debug.debug_file) << leaf -> getFullString() << "\n";
			}

			if(extensions.size() > 0)
			{
				updateLeaves(newLeaves, extensions, *iter,currLeavesNum);
				break;
			}
			isReduced = false;
			m_min_SA_threshold--;
			count++;
		}
		m_min_SA_threshold += count;

		if (minTotalcount >= totalcount)
		{
			minTotalcount = totalcount;
		}

		++iter;
		++currLeavesNum;
	}
}

// Update the leaves after the extension is successful.
void LongReadSelfCorrectByOverlap::updateLeaves(leafList& newLeaves,extArray& extensions,leafInfo& leaf,size_t currLeavesNum)
{
	SAIOverlapNode3* pNode    = leaf.leafNodePtr;
	if(extensions.size() == 1)
	{
		// Single extension, do not branch
			pNode->extend(extensions.front().SearchLetters);
			newLeaves.emplace_back(pNode, leaf, extensions.front(), currLeavesNum);
	}
	else if(extensions.size() > 1)
	{
		// Branch
		for(size_t i = 0; i < extensions.size(); ++i)
		{
			SAIOverlapNode3* pChildNode = pNode->createChild(extensions[i].SearchLetters);
			//inherit accumulated kmerCount from parent
				pChildNode->addKmerCount( pNode->getKmerCount() );
/*
			// Remove the bug while changing the seed by FM-Extend
				if (i != 0)
					(pChildNode -> resultindex).first = -1;
//*/
			newLeaves.emplace_back(pChildNode, leaf, extensions[i], currLeavesNum);
		}
	}
}

// Compute the error rates and keep the leaves whose error rates <= expected error rate.
bool LongReadSelfCorrectByOverlap::PrunedBySeedSupport(leafList& newLeaves)
{
	// the seed index in m_TerminatedIntervals for m_currentLength
	// the m_currentLength is the same for all leaves
	// which is used as the central index within the m_maxIndelSize window
	size_t currSeedIdx = m_currentLength-m_seedSize;

	//        ---	seed size =3. seed dist = 1;
	//         --*
	//          -*-
	//           *--
	//            ---
	// *: SNP or indel errors, ---: seed size
	// Erase the leaf if no feasible seeds are found within seedSize+maxIndelSize.
	size_t indelOffset = m_seedSize+m_maxIndelSize;

	// Compute the range of small and large indices for tolerating m_maxIndelSize
	size_t smallSeedIdx = currSeedIdx <= indelOffset ? 0 : currSeedIdx - indelOffset;
	size_t largeSeedIdx = (currSeedIdx+indelOffset) >= (m_query.length()-m_seedSize)?
						  (m_query.length()-m_seedSize):currSeedIdx+indelOffset;

	// check range of last seed and find new seeds for each interval
	leafList::iterator iter = newLeaves.begin();
	while(iter != newLeaves.end())
	{
		bool isNewSeedFound = false;
		SAIOverlapNode3* leaf = (*iter).leafNodePtr;

		if (m_currentLength - leaf -> lastOverlapLen > m_seedSize
		||  m_currentLength - leaf -> lastOverlapLen <= 1 )
		{
			size_t preSeedIdx = leaf -> lastSeedIdx;
			// search for matched new seeds
			isNewSeedFound = isSupportedByNewSeed(leaf, smallSeedIdx, largeSeedIdx);

			// lastSeedIdxOffset records the offset between lastSeedIdx and currSeedIdx when first match is found
			if(isNewSeedFound)
			{
				if( currSeedIdx + leaf -> lastSeedIdxOffset - preSeedIdx > m_seedSize )
					leaf -> numRedeemSeed += (m_seedSize-1)*m_PacBioErrorRate;

				leaf -> lastSeedIdxOffset = (int) leaf->lastSeedIdx - (int)currSeedIdx;
			}
			else
			{
				// If the seed extension is stopped by SNP or indel error for the 1st time
				// increment the error number in order to distinguish two separate seeds
				// and one larger consecutive seed during error rate computation
				if     ( (currSeedIdx + leaf->lastSeedIdxOffset - leaf->lastSeedIdx) % m_seedSize == 1 )
					leaf->numOfErrors ++;
				else if( (currSeedIdx + leaf->lastSeedIdxOffset - leaf->lastSeedIdx) > m_seedSize  - 1 )
					leaf->numRedeemSeed += 1-m_PacBioErrorRate;
			}
		}

		else
			leaf->numRedeemSeed += 1-m_PacBioErrorRate;

		double currErrorRate = computeErrorRate(leaf);

		// speedup by skipping dissimilar reads
		// This is the 2nd filter less reliable than the 1st one
		if((!m_Debug.ref.isSpecifiedPath) && (currErrorRate > m_errorRate)) //testcw
		{
			iter = newLeaves.erase(iter);
			continue;
		}

		iter++;
	}

	return true;
}

// Find the matched k-mer with the current extension sequence in the raw read.
bool LongReadSelfCorrectByOverlap::isSupportedByNewSeed(SAIOverlapNode3* currNode, size_t smallSeedIdx, size_t largeSeedIdx)
{
	// If there is mismatch/indel, jump to the next m_seedSize/m_seedDist, and 1 otherwise.
	size_t seedIdxOffset = currNode->lastOverlapLen < m_currentLength-m_seedSize?
							m_seedSize:m_currentLength - currNode->lastOverlapLen;

	// search for new seed starting from last matched seed or smallSeedIdx
	size_t startSeedIdx = std::max(smallSeedIdx, currNode->lastSeedIdx+seedIdxOffset);

	bool isNewSeedFound = false;
	BWTInterval currFwdInterval = currNode->fwdInterval;
	BWTInterval currRvcInterval = currNode->rvcInterval;

	// Binary search for new seeds using Query interval tree
	std::vector<TreeInterval<size_t> > resultsFwd, resultsRvc;
	if(currFwdInterval.isValid())
		m_fwdIntervalTree.findOverlapping(currFwdInterval.lower, currFwdInterval.upper, resultsFwd);
	if(currRvcInterval.isValid())
		m_rvcIntervalTree.findOverlapping(currRvcInterval.lower, currRvcInterval.upper, resultsRvc);
	int minIdxDiff = 10000;
	size_t currSeedIdx = m_currentLength-m_seedSize;
	for(size_t i=0 ; i<resultsFwd.size() || i<resultsRvc.size() ; i++)
	{
		if( currFwdInterval.isValid() &&
			i<resultsFwd.size() &&
			resultsFwd.at(i).value >= startSeedIdx &&
			resultsFwd.at(i).value <= largeSeedIdx )
		{
			// update currNode members
			if(std::abs((int)resultsFwd.at(i).value - (int)currSeedIdx) < minIdxDiff)
			{
				currNode->lastSeedIdx = resultsFwd.at(i).value;

				// query overlap may shift due to indels
				currNode->queryOverlapLen = resultsFwd.at(i).value+m_seedSize;
				minIdxDiff = std::abs((int)resultsFwd.at(i).value - (int)currSeedIdx);
			}
			// lastOverlapLen records the overlap length of last hit
			currNode->lastOverlapLen = m_currentLength;
			// currOverlapLen is always identical to m_currentLength
			currNode->currOverlapLen = m_currentLength;
			isNewSeedFound = true;
		}
		else if( currRvcInterval.isValid() &&
			i<resultsRvc.size() &&
			resultsRvc.at(i).value >= startSeedIdx &&
			resultsRvc.at(i).value <= largeSeedIdx )
		{

			// update currNode members
			if(std::abs( (int)currSeedIdx - (int)resultsRvc.at(i).value ) < minIdxDiff)
			{
				currNode->lastSeedIdx = resultsRvc.at(i).value;

				// query overlap may shift due to indels
				currNode->queryOverlapLen = resultsRvc.at(i).value+m_seedSize;
				minIdxDiff = std::abs((int)resultsRvc.at(i).value - (int)currSeedIdx);
			}
			// lastOverlapLen records the overlap length of last hit
			currNode->lastOverlapLen = m_currentLength;
			// currOverlapLen is always identical to m_currentLength
			currNode->currOverlapLen = m_currentLength;
			isNewSeedFound = true;
		}
	}

	if(isNewSeedFound)
		currNode->totalSeeds++;
	return isNewSeedFound;
}

// Compute the error rate in the leaf.
double LongReadSelfCorrectByOverlap::computeErrorRate(SAIOverlapNode3* currNode)
{
	// Compute accuracy via matched length in both query and subject

	double matchedLen = (double)currNode->totalSeeds + m_seedSize-1;

	// SNP and indel over-estimate the unmatched lengths across error, ---*---
	// Restore the unmatched region via numOfErrors, which is still over-estimated

	matchedLen += currNode->numRedeemSeed;

	double totalLen = (double)currNode->currOverlapLen ;

	double unmatchedLen = totalLen - matchedLen;

	double  currErrorRate =  unmatchedLen/totalLen;
	currNode->GlobalErrorRateRecord.push_back(currErrorRate);

	if(currNode->GlobalErrorRateRecord.size() >= m_localSimilarlykmerSize)
	{
		size_t totalsize = currNode->GlobalErrorRateRecord.size();

		currErrorRate = ( currErrorRate*totalLen-currNode->GlobalErrorRateRecord.at( totalsize - m_localSimilarlykmerSize)*(totalLen - m_localSimilarlykmerSize) )/m_localSimilarlykmerSize;
	}
	currNode->LocalErrorRateRecord.push_back(currErrorRate);
	return currErrorRate;
}

// Attempt to the extension in the current leaf and get the extension information.
extArray LongReadSelfCorrectByOverlap::getFMIndexExtensions(const leafInfo& currLeaf,const bool printDebugInfo)
{
	SAIOverlapNode3* leaf = currLeaf.leafNodePtr;
	extArray output;
		output.reserve(4);
	extArray totalExt;
		totalExt.reserve(4);

	size_t IntervalSizeCutoff = m_min_SA_threshold;    //min freq at fwd and rvc bwt, >=3 is equal to >=2 kmer freq

	totalcount = 0;
	int maxfreqsofleaf = 0;
	dominantBase maxFreqOfLeaves(currLeaf.tailLetter);
/*
	if(m_Debug.isDebug)
		std::cout   << leaf->getFullString() <<" || Local Error Rate: "
					<< leaf->LocalErrorRateRecord.back()  <<"\n"
					<< leaf->getSuffix(m_currentKmerSize) << " || k-mer freqs: "
					<< currLeaf.kmerFrequency << std::endl;
*/
	for(int i = 1; i < BWT_ALPHABET::size; ++i) //i=A,C,G,T
	{
		//update forward Interval using extension b
			char b = BWT_ALPHABET::getChar(i);
			BWTInterval fwdProbe = leaf->fwdInterval;
			if(fwdProbe.isValid())
				BWTAlgorithms::updateInterval(fwdProbe,b,m_pRBWT);

		//update reverse complement Interval using extension rcb
			char rcb=BWT_ALPHABET::getChar(5-i);
			BWTInterval rvcProbe = leaf->rvcInterval;
			if(rvcProbe.isValid())
				BWTAlgorithms::updateInterval(rvcProbe,rcb,m_pBWT);

		FMidx currExt = FMidx(b,fwdProbe,rvcProbe);

		totalcount += currExt.getKmerFrequency();

		if(m_Debug.isDebug && printDebugInfo)
		{
			*(m_Debug.debug_file) << currExt.getKmerFrequency();
			if(i < BWT_ALPHABET::size-1)
				*(m_Debug.debug_file) << "|";
			else
				*(m_Debug.debug_file) << "&";
		}

		maxFreqOfLeaves.setFreq(currExt);

		totalExt.push_back(currExt);

	}// end of ACGT

	m_maxfreqs = std::max(m_maxfreqs,totalcount);
	maxfreqsofleaf = maxFreqOfLeaves.getMaxFreq();

	for(int i = 1; i < BWT_ALPHABET::size; ++i)
	{
		// Get the information
		size_t kmerFreq         = totalExt.at(i-1).getKmerFrequency();
		BWTInterval fwdInterval = totalExt.at(i-1).getFwdInterval();
		BWTInterval rvcInterval = totalExt.at(i-1).getRvcInterval();

		// Compute the k-mer ratio arguments
		const double kmerRatioNotPass = 2;
		double kmerRatioCutoff = 0;
		double kmerRatio = (double) kmerFreq/(double)maxfreqsofleaf;

		char b = BWT_ALPHABET::getChar(i);

		bool isHomopolymer = (currLeaf.tailLetterCount >= 3);
		bool isMatchedBy5mer = ismatchedbykmer(fwdInterval,rvcInterval);

		bool isDominant     = maxFreqOfLeaves.isDominant();
		bool isFreqPass     = kmerFreq   >= IntervalSizeCutoff;
		bool isLowCoverage  = totalcount >= IntervalSizeCutoff+2;
		bool isRepeat       = maxfreqsofleaf > 100;
		bool isHighlyRepeat = maxfreqsofleaf > 150;
		bool isLowlyRepeat  = maxfreqsofleaf >  50;

		// matched case
			if  ( isMatchedBy5mer &&  isHighlyRepeat )
				kmerRatioCutoff = 0.125;
		else if ( isMatchedBy5mer &&  isLowlyRepeat  )
				kmerRatioCutoff = 0.2  ;
		// unmatched case
		else if ( isFreqPass )
				kmerRatioCutoff = 0.25 ;
		else if ( isLowCoverage )
				kmerRatioCutoff = 0.6  ;
		else
				kmerRatioCutoff = kmerRatioNotPass;

		// Homopolymer case
			if  ( isHomopolymer   &&  isRepeat )
				kmerRatioCutoff = std::max(kmerRatioCutoff,0.3);
		else if ( isHomopolymer )
				kmerRatioCutoff = std::max(kmerRatioCutoff,0.6);
/*
			if (isDominant)
				kmerRatioCutoff = std::max(kmerRatioCutoff,0.4);
//*/
		if(m_Debug.isDebug && printDebugInfo)
		{
			*(m_Debug.debug_file) << kmerRatio;
			if(i < BWT_ALPHABET::size-1)
				*(m_Debug.debug_file) << "|";
			else
				*(m_Debug.debug_file) << "&";
		}

		if( kmerRatio >=   kmerRatioCutoff )
		{
			// extend to b
				output.emplace_back(b,fwdInterval,rvcInterval);
		}
	}
	if(m_Debug.ref.isSpecifiedPath)
	{
		if (!(m_step_number == 1 || m_Debug.ref.isMatchBase(currLeaf.tailLetter, m_step_number-2)))
		{
			output.clear();
		}
		else if(!output.empty() || !printDebugInfo)
		{
			// output is empty, and the k-mer size reduced
			// output is not empty

			output.clear();
			for(int i = 1; i < BWT_ALPHABET::size; ++i)
			{
				BWTInterval fwdInterval = totalExt.at(i-1).getFwdInterval();
				BWTInterval rvcInterval = totalExt.at(i-1).getRvcInterval();
				char b = BWT_ALPHABET::getChar(i);
				output.emplace_back(b,fwdInterval,rvcInterval);
			}
		}
	}

	if(m_Debug.isDebug && printDebugInfo)
	{
		*(m_Debug.debug_file)
			<< maxfreqsofleaf  << "|"
			<< totalcount      << "&"
			<< m_extSeeds.source.isRepeat << "|"
			<< m_extSeeds.target.isRepeat << "\n";
	}

	return output;
}

// Determine if the current sequence is matched 5-mer raw read.
bool LongReadSelfCorrectByOverlap::ismatchedbykmer(BWTInterval currFwdInterval,BWTInterval currRvcInterval)
{
	bool match = false;

	// Binary search for new seeds using Query interval tree
	std::vector<TreeInterval<size_t> > resultsFwd, resultsRvc;
	if(currFwdInterval.isValid())
		m_fwdIntervalTree2.findOverlapping(currFwdInterval.lower, currFwdInterval.upper, resultsFwd);
	if(currRvcInterval.isValid())
		m_rvcIntervalTree2.findOverlapping(currRvcInterval.lower, currRvcInterval.upper, resultsRvc);
	size_t startSeedIdx = std::max((int)m_currentLength - (int)m_maxIndelSize,0);
	size_t largeSeedIdx = m_currentLength + m_maxIndelSize;

	for(size_t i=0 ; i<resultsFwd.size() || i<resultsRvc.size() ; i++)
	{
		if( currFwdInterval.isValid() &&
			i<resultsFwd.size() &&
			resultsFwd.at(i).value >= startSeedIdx &&
			resultsFwd.at(i).value <= largeSeedIdx )
		{
			match = true;
			break;
		}
		else if( currRvcInterval.isValid() &&
			i<resultsRvc.size() &&
			resultsRvc.at(i).value >= startSeedIdx &&
			resultsRvc.at(i).value <= largeSeedIdx )
		{
			match = true;
			break;
		}
	}

	return match;
}

// Check for leaves whose extension has terminated.
// If the leaf has terminated, the walked string and coverage is pushed to the result vector
bool LongReadSelfCorrectByOverlap::isTerminated(SAIntervalNodeResultVector& results)
{
	bool found = false;

	for(leafList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
	{
		SAIOverlapNode3* leaf = (*iter).leafNodePtr;
		BWTInterval currfwd = leaf -> fwdInterval;
		BWTInterval currrvc = leaf -> rvcInterval;

		if (!m_Debug.ref.isSpecifiedPath)
			assert(currfwd.isValid() || currrvc.isValid());

		//The current SA interval stands for a string >= terminating kmer
		//If terminating kmer is a substr, the current SA interval is a sub-interval of the terminating interval
		bool isFwdTerminated = false;
		bool isRvcTerminated = false;
		for(size_t i = std::max(leaf -> resultindex.second,0);i < m_fwdTerminatedInterval.size(); i++)
		{

			isFwdTerminated=currfwd.isValid() && currfwd.lower >= m_fwdTerminatedInterval.at(i).lower
							&& currfwd.upper <= m_fwdTerminatedInterval.at(i).upper;
			isRvcTerminated=currrvc.isValid() && currrvc.lower >= m_rvcTerminatedInterval.at(i).lower
							&& currrvc.upper <= m_rvcTerminatedInterval.at(i).upper;

			if(isFwdTerminated || isRvcTerminated)
			{
				std::string STNodeStr = leaf->getFullString();
				if (m_targetSeed.length() > m_minOverlap)
					STNodeStr += m_targetSeed.substr(i+m_minOverlap);

				SAIntervalNodeResult STresult;
				STresult.thread=STNodeStr;
				STresult.SAICoverage = leaf->getKmerCount();
				STresult.errorRate   = leaf->GlobalErrorRateRecord.back();
				STresult.SAIntervalSize = (currfwd.upper-currfwd.lower+1);

				if( leaf->resultindex.first == -1 )
				{
					results.push_back(STresult);
					leaf->resultindex = std::make_pair(results.size(),i);
				}
				else
				{
					results.at( leaf->resultindex.first-1 ) = STresult;
					leaf->resultindex = std::make_pair(leaf->resultindex.first,i);
				}

				found =  true;
			}
		}

	}

	return found;
}

void LongReadSelfCorrectByOverlap::printErrorRate(leafList& currLeaves)
{
	std::cout << std::fixed << std::setprecision(5);

	for(auto& iter : currLeaves)
	{
		SAIOverlapNode3* leaf = iter.leafNodePtr;
		std::cout << 100*(leaf->LocalErrorRateRecord.back()) << "\t";
	}

	std::cout << std::endl;
}