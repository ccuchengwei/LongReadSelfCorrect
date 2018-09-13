///----------------------------------------------
// Copyright 2016 National Chung Cheng University
// Written by Yao-Ting Huang & Ping-Yeh Chen
// Released under the GPL
//-----------------------------------------------
//
// LongReadSelfCorrectByFMExtend - A fast overlap filter using locality-sensitive backward search.
//
//
#include <iomanip>
#include "LongReadCorrectByFMExtend.h"
#include "KmerFeature.h"
#include "BWTAlgorithms.h"
#include "stdaln.h"
#include "LongReadOverlap.h"

thread_local std::vector<double> LongReadSelfCorrectByFMExtend::freqsOfKmerSize;

// Class: SAIOverlapTree
LongReadSelfCorrectByFMExtend::LongReadSelfCorrectByFMExtend
				(
					const size_t caseNumber,
					const seedPair& extSeeds,
					const std::string& strBetweenSrcTarget,
					int disBetweenSrcTarget,
					size_t initkmersize,
					size_t maxResetSize,
					const FMextendParameters params,
					size_t min_SA_threshold,
					const debugExtInfo debug,
					double errorRate,
					size_t repeatFreq,
					size_t localSimilarlykmerSize
				):
					m_case_number(caseNumber),
					m_extSeeds(extSeeds),
					m_sourceSeed(extSeeds.source.seq),
					m_strBetweenSrcTarget(strBetweenSrcTarget),
					m_targetSeed(extSeeds.target.seq),
					m_disBetweenSrcTarget(disBetweenSrcTarget),
					m_initkmersize(initkmersize),
					m_minOverlap(params.minKmerLength),
					m_maxResetSize(maxResetSize),
					m_pBWT(params.indices.pBWT),
					m_pRBWT(params.indices.pRBWT),
					m_PBcoverage(params.PBcoverage),
					m_min_SA_threshold(min_SA_threshold),
					m_errorRate(errorRate),
					m_maxLeaves(params.maxLeaves),
					m_matchedSize(params.idmerLength),
					m_repeatFreq(repeatFreq),
					m_localSimilarlykmerSize(localSimilarlykmerSize),
					m_PacBioErrorRate(params.ErrorRate),
					m_Debug(debug),
					fieldNum(0),
					m_currLeavesTotalSize(1),
					m_currLeavesKeptSize (1),
					m_nextLeavesTotalSize(0),
					m_nextLeavesKeptSize (0),
					isSeparatedInitKmer(false),
					m_totalExtNum(1),
					m_accumLeaves(1)
{
	m_extSeeds.reduceSourceBy(m_sourceSeed.length()-m_initkmersize);

	//if distance < 100 ,use const indel size
		if (m_disBetweenSrcTarget > 100)
			m_maxIndelSize  =  m_disBetweenSrcTarget * 0.2;
		else
			m_maxIndelSize  =  20;

	//frequencies of correspond k
		// The array don't change at everytime FM-extension.
		// Thus, it only needed declared at first FM-extension.
		if (m_case_number == 1)
		{
			const int maxIndex = 100 + 1;
			LongReadSelfCorrectByFMExtend::freqsOfKmerSize.assign(maxIndex,0);
			for(int i = m_minOverlap ; i < maxIndex ; i++)
				LongReadSelfCorrectByFMExtend::freqsOfKmerSize.at(i)
											= pow(1 - m_PacBioErrorRate, i) * m_PBcoverage;
		}

	// PacBio reads are longer than real length due to insertions
		m_maxLength =          (1.2*(m_disBetweenSrcTarget+10))+2*m_initkmersize     ;
		m_minLength = std::max((0.8*(m_disBetweenSrcTarget-20))+2*m_initkmersize,0.0);

	//initialRootNode
		initialRootNode(m_extSeeds.source.seq);

	// push new node into roots and leaves vector
		m_RootNodes .push_back   (m_pRootNode);
		m_currLeaves.emplace_back(m_pRootNode);

		totalPathInfo.reserve(m_maxLength+2);
		setPathInfo(totalPathInfo, m_pRootNode, m_extSeeds.source.seq, isSeparatedInitKmer);

		if (isSeparatedInitKmer)
			m_pathOffset = m_initkmersize - 1;
		else
			m_pathOffset = 0;

//*
	if(m_Debug.isDebug)
	{
		(*m_Debug.debug_finalSeq)
					<< m_case_number << "\tBE: "
					<< m_extSeeds.source.seq << " "
					<< m_extSeeds.target.seq << " "
					<< m_pRootNode->fwdInterval.size() + m_pRootNode->rvcInterval.size() << "|"
					<< disBetweenSrcTarget   << " ("
					<< m_minLength           << ", "
					<< m_maxLength           << ") "
					<< m_maxIndelSize        << "\n";
	}
//*/

	// initialize the ending SA intervals with kmer length = m_minOverlap
		// const size_t tailIndex = 0; // Don't change the seed
		const size_t tailIndex = m_targetSeed.length()-m_minOverlap;
		for(size_t i =0 ;i <= tailIndex; i++)
		{
			std::string endingkmer = m_targetSeed.substr(i, m_minOverlap);

			m_fwdTerminatedInterval.push_back(BWTAlgorithms::findInterval(m_pRBWT, reverse(endingkmer)));
			m_rvcTerminatedInterval.push_back(BWTAlgorithms::findInterval(m_pBWT, reverseComplement(endingkmer)));
		}

    // build overlap tree
		m_query = m_extSeeds.source.seq + m_strBetweenSrcTarget + m_targetSeed;
		// build overlap tree to determine the error rate
			buildOverlapbyFMindex(m_fwdIntervalTree ,m_rvcIntervalTree ,m_matchedSize);
		// build overlap tree to match 5-mer
			buildOverlapbyFMindex(m_fwdIntervalTree2,m_rvcIntervalTree2,5);

}

LongReadSelfCorrectByFMExtend::~LongReadSelfCorrectByFMExtend()
{
	for (auto iter = m_RootNodes.begin(); iter != m_RootNodes.end(); ++iter)
		delete *iter;

	m_RootNodes.clear();
}

// Initialize the root node
void LongReadSelfCorrectByFMExtend::initialRootNode(const std::string& beginningkmer)
{
	// create one root node
		m_pRootNode = new SAIOverlapNode3(&m_sourceSeed, NULL);
		m_pRootNode -> setRoot(beginningkmer, m_initkmersize, m_matchedSize, m_pBWT, m_pRBWT);

		m_maxfreqs = m_pRootNode->fwdInterval.size() + m_pRootNode->rvcInterval.size();
		m_currentLength = m_currentKmerSize = m_initkmersize;
}

//Build the overlap tree
void LongReadSelfCorrectByFMExtend::buildOverlapbyFMindex(IntervalTree<size_t>& fwdIntervalTree,IntervalTree<size_t>& rvcIntervalTree,const int& overlapSize)
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
int LongReadSelfCorrectByFMExtend::extendOverlap(FMWalkResult2& FMWResult)
{
	SAIntervalNodeResultVector results;
	m_step_number = 1;

	//Overlap extension via FM-index walk
	while((m_currLeavesKeptSize != 0) && (m_currLeavesKeptSize <= m_maxLeaves) && (m_currentLength <= m_maxLength))
	{
		// ACGT-extend the leaf nodes via updating existing SA interval
			m_nextLeaves.clear();
			extendLeaves();

		// Use overlap tree to trim branch
			PrunedByAbsErrorRate();

		// Change old leaves to the new leaves
			m_currLeaves.clear();
			m_currLeaves          = m_nextLeaves;
			m_currLeavesKeptSize  = m_nextLeavesKeptSize ;
			m_currLeavesTotalSize = m_nextLeavesTotalSize;

		// If the length is enough,
		// then confirm whether the candidate path reaches the target seed.
			if(m_currentLength >= m_minLength)
				isTerminated(results);

		// Change to the next step
			m_step_number++;

		// The debug mode to force to extend into the correct path.
			if (m_Debug.ref.isSpecifiedPath)
			{
				if (m_Debug.ref.isIllegalLen(m_step_number-1))
				{
					m_currLeaves.clear();
					m_currLeavesKeptSize  = m_currLeavesTotalSize = 0;
					break;
				}
			}
	}

	// reach the terminal k-mer
	if(results.size() > 0)
		return findTheBestPath(results, FMWResult); // find the path with maximum match percent or kmer coverage

	// Did not reach the terminal kmer
	if      (m_currLeavesKeptSize == 0          ) // high error
		return -1;
	else if (m_currentLength      >  m_maxLength) // exceed search depth(reach to the error path)
		return -2;
	else if (m_currLeavesKeptSize >  m_maxLeaves) // too much repeats
		return -3;
	else
		return -4;
}

// Find the path where is the min error rate.
int LongReadSelfCorrectByFMExtend::findTheBestPath(const SAIntervalNodeResultVector& results, FMWalkResult2& FMWResult)
{
	double minErrorRate = 1;
/*
	// Declare variables to prepare for printing debug v.2.
		size_t realCandiNum = 1;
		size_t mismatchSNum = 1;
//*/
	for (size_t i = 0 ; i < results.size() ;i++)
	{
		const std::string& candidateSeq = results.at(i).thread;
//*
		//  Print debug v.1.
			if(m_Debug.isDebug)
			{
				(*m_Debug.debug_finalSeq)
								<< "Final Seqs: "   << results.at(i).thread
								<< "\tError Rate: " << results.at(i).errorRate << std::endl;
			}
//*/
		if( results.at(i).errorRate < minErrorRate )
		{
/*
			// Initialize variables to prepare for printing debug v.2.
				realCandiNum = 1;
				mismatchSNum = 1;
//*/
			// Store the results
				minErrorRate = results.at(i).errorRate;
				FMWResult.mergedSeq = candidateSeq;
				minTotalcount = results.at(i).SAIntervalSize;
		}
/*
		else if (results.at(i).errorRate == minErrorRate )
		{
			// Compute variables to prepare for printing debug v.2.
				realCandiNum++;
				if (m_Debug.isDebug)
				{
					if (candidateSeq.length() == FMWResult.mergedSeq.length())
					{
						mismatchSNum++;
					}
				}
		}
//*/
	}
/*
	//  Print debug v.2.
		if (m_Debug.isDebug)
		{
			(*m_Debug.debug_finalSeq)
							<< m_Debug.readID        << "\t"
							<< m_case_number         << "\t"
							<< FMWResult.mergedSeq.length() << "\t"
							<< realCandiNum          << "\t"
							<< mismatchSNum          << "\t"
							<< m_accumLeaves         << std::endl;
		}
//*/
	if(FMWResult.mergedSeq.length() != 0)
		return 1;
	return -4;
}

// Extend all Leaves by FM-index
void LongReadSelfCorrectByFMExtend::extendLeaves()
{
	// Set the debug parameters
		// the first  bit: is print the info. (original or reducedSize)
		// the second bit: reduce the threshold by "extendLeaves".
		// the third  bit: reduce the threshold by "attempToExtend".
			debugBits_t debugBits(0); // 0 = (000)2
			std::string debugStr("");
			std::string debugStrForSizeReduced("");


	//resize if length too long
		if(m_currentKmerSize > m_maxResetSize)
			refineSAInterval(m_currLeaves, m_maxResetSize);

	// Set the leaves information
		const   leafTable_t& lastPathInfo = totalPathInfo.back();
				leafTable_t  currPathInfo;
				// the max branches = 4 per leaf.
					currPathInfo.reserve(4 * m_currLeavesKeptSize);
					m_nextLeaves.reserve(4 * m_currLeavesKeptSize);

	// The first attempt to extend the leaves
		PrunedByRelErrorRate(lastPathInfo);
		attempToExtend(debugBits, debugStr, lastPathInfo, currPathInfo);

	if( m_nextLeavesKeptSize == 0 ) //level 1 reduce size
	{
		size_t LowerBound = std::max(m_currentKmerSize - 2, m_minOverlap);
		size_t ReduceSize = SelectFreqsOfrange(LowerBound, m_currentKmerSize, m_currLeaves);

		// Set the debug information
			debugBits.set(0, (m_currentKmerSize != ReduceSize));

		refineSAInterval(m_currLeaves, ReduceSize);
		attempToExtend(debugBits , debugStrForSizeReduced, lastPathInfo, currPathInfo);

		if( m_nextLeavesKeptSize == 0 )//level 2 reduce threshold
		{
			// Set the debug information
				if (m_Debug.isDebug)
				{
					debugStrForSizeReduced.clear();
					debugBits.set(1,true);
				}
			m_min_SA_threshold--;
			attempToExtend(debugBits , debugStrForSizeReduced, lastPathInfo, currPathInfo);
			m_min_SA_threshold++;
		}

		// Set the debug information
			if (m_Debug.isDebug)
			{
				if (debugBits[0])
					debugStr += debugStrForSizeReduced;
				else
					debugStr  = debugStrForSizeReduced;
			}
	}

	// Set the debug information
		if (m_Debug.isDebug)
			*(m_Debug.debug_file) << debugStr;

	m_totalExtNum += m_currTotalExtNum;

	// Extension succeed
		if( m_nextLeavesKeptSize != 0 )
		{
			m_currentLength++;
			m_currentKmerSize++;
			totalPathInfo.emplace_back(currPathInfo);

			// if frequency are low , relax it to prepare for the next extension.
				if( isInsufficientFreqs() )
				{
					size_t LowerBound = std::max(m_currentKmerSize - 2, m_minOverlap);
					size_t ReduceSize = SelectFreqsOfrange(LowerBound,m_currentKmerSize,m_nextLeaves);
					refineSAInterval(m_nextLeaves,ReduceSize);
				}

		}

}

// Determine the size of the reduced k-mer in leaves when the max freq of them approach to the expected one.
size_t LongReadSelfCorrectByFMExtend::SelectFreqsOfrange(const size_t LowerBound, const size_t UpperBound, LeavesArray_t& leaves)
{
	extArray maxKmerArray;
	int tempmaxfmfreqs = 0;

	for(auto& leaf : leaves)
	{
		if ( leaf -> isRemoved() )
			continue;

		std::string    maxKmer = leaf -> getSuffix(UpperBound);

		std::string startkmer  = maxKmer.substr(UpperBound - LowerBound); //  string of lower bound kmer size

		BWTInterval Fwdinterval = BWTAlgorithms::findInterval(m_pBWT, startkmer);
		BWTInterval Rvcinterval = BWTAlgorithms::findInterval(m_pRBWT, reverseComplement(reverse(startkmer)));

		maxKmerArray.emplace_back(maxKmer,Fwdinterval,Rvcinterval);
		FMidx_t& currKmer = maxKmerArray.back();

		if(currKmer.getKmerFrequency() > tempmaxfmfreqs ) //check interval size
			tempmaxfmfreqs = currKmer.getKmerFrequency();

	}

	if( tempmaxfmfreqs - (int)LongReadSelfCorrectByFMExtend::freqsOfKmerSize.at(LowerBound) < 5 ) return LowerBound;

	for(size_t i=1 ; i <= UpperBound - LowerBound; i++ )
	{
		tempmaxfmfreqs = 0;
		for(size_t j = 0; j < maxKmerArray.size(); j++)
		{
			std::string startkmer   = maxKmerArray.at(j).SearchLetters.substr(UpperBound - LowerBound - i);
			BWTInterval Fwdinterval = maxKmerArray.at(j).getFwdInterval();
			BWTInterval Rvcinterval = maxKmerArray.at(j).getRvcInterval();

			char b   = startkmer.front();//b:base
			char rcb = complement(b);
			BWTAlgorithms::updateInterval(Fwdinterval,  b,m_pBWT );
			BWTAlgorithms::updateInterval(Rvcinterval,rcb,m_pRBWT);

			maxKmerArray.at(j).setInterval(Fwdinterval,Rvcinterval);

			if(maxKmerArray.at(j).getKmerFrequency() > tempmaxfmfreqs) //check interval size
				tempmaxfmfreqs = maxKmerArray.at(j).getKmerFrequency();

		}

		if( tempmaxfmfreqs - (int)LongReadSelfCorrectByFMExtend::freqsOfKmerSize.at(LowerBound + i) < 5 ) return LowerBound + i ;
	}

	return UpperBound;
}

// Determine if frequencies in leaves are too low to perform the next extension.
bool LongReadSelfCorrectByFMExtend::isInsufficientFreqs()
{
	size_t highfreqscount = 0;
	for(auto& leaf : m_nextLeaves)
	{
		size_t highfreqThreshold = m_PBcoverage > 60 ? (size_t)(m_PBcoverage/60)*3 : 3;

		if((leaf -> isKept()) && (leaf -> getLastKmerCount()) > highfreqThreshold)
			highfreqscount++;
	}

	if( highfreqscount == 0 )
		return true;
	else if ( highfreqscount <= 2 && m_nextLeavesKeptSize >= 5 )
		return true;
	else if ( highfreqscount <= 1 && m_nextLeavesKeptSize >= 3 )
		return true;
	return false;
}

// Refine SA intervals of each leave with a new k-mer
void LongReadSelfCorrectByFMExtend::refineSAInterval(LeavesArray_t& leaves, const size_t newKmerSize)
{
	for(auto& leaf : leaves)
	{
		// Reset the SA intervals using newKmerSize
		if (leaf -> isRemoved())
			continue;

			std::string reducedKmer = leaf -> getSuffix(newKmerSize);

			leaf -> fwdInterval = BWTAlgorithms::findInterval(m_pRBWT,           reverse(reducedKmer));
			leaf -> rvcInterval = BWTAlgorithms::findInterval(m_pBWT , reverseComplement(reducedKmer));
	}

	m_currentKmerSize = newKmerSize;
}

// Keep the leaves whose highly-correct rates are relative to the others.
void LongReadSelfCorrectByFMExtend::PrunedByRelErrorRate(const leafTable_t& lastPathInfo)
{
	// Compute the min error rate
		m_minMatchedErrorRate = 1.0;
		for(auto& leaf : m_currLeaves)
		{
			if((leaf -> isKept()) && lastPathInfo.at(leaf).localErrorRate < m_minMatchedErrorRate)
				m_minMatchedErrorRate = lastPathInfo.at(leaf).localErrorRate;
		}

	// Compute the errorRateDiff to trim leaves whose error rates relative to the others is high.
		for(auto& leaf : m_currLeaves)
		{
			if (m_Debug.ref.isSpecifiedPath)
				break;

			double errorRateDiff  = lastPathInfo.at(leaf).localErrorRate - m_minMatchedErrorRate;
			if  (leaf -> isKept())
			{
				if (
							(errorRateDiff > 0.05 && m_currentLength > m_localSimilarlykmerSize/2)
						||  (errorRateDiff > 0.1  && m_currentLength > 15)
					)
					{
						leaf -> removeType = -1;
						m_currLeavesKeptSize--;
					}
			}
		}
}

// And attempt to extend those leaves.
void LongReadSelfCorrectByFMExtend::attempToExtend(const debugBits_t& debugBits, std::string& debugStr, const leafTable_t& lastPathInfo, leafTable_t& currPathInfo)
{
	// Set the debug parameters
		debugBits_t  pDebugBits (debugBits);
		std::string  debugSubstr("");
		std::stringstream strBuffer;

		m_maxfreqs = 0;
		minTotalcount = 10000000;
		m_currTotalExtNum  = 0;

		double maximumAverfreqs = 0.0;
		size_t absolutePruning  = 0;
		size_t relativePruning  = 0;
		if (m_Debug.isDebug)
		{
			for(auto& leaf : m_currLeaves)
			{
				kmerFreq_t allKmer_freq = leaf -> getKmerCount();
				double     aveKmer_freq = (double) allKmer_freq/(double) m_currentLength;

				if (leaf -> isKept())
				{
					if( maximumAverfreqs < aveKmer_freq )
						maximumAverfreqs = aveKmer_freq ;
				}
				else if (leaf -> removeType ==  1)
					absolutePruning++;
				else if (leaf -> removeType == -1)
					relativePruning++;
			}
		}

	size_t currLeavesNum = 0;
	for (auto leaf : m_currLeaves)
	{
		extArray extensions;
		int count = 0;
		pDebugBits.set(0, debugBits[0]);
		while(count < 2)
		{
			if	(
					count == 1
					&&
					(
							leaf -> isRemoved()
						||  lastPathInfo.at(leaf).localErrorRate != m_minMatchedErrorRate
						|| m_currLeavesKeptSize == 1
					)
				)
				break; // Relax the leaf with the minimum error rate.
			else
			{
				// Reset the debug strings
					strBuffer.str("");
					strBuffer.clear();
			}

			// Set the debug strings
				if (m_Debug.isDebug)
				{
					pDebugBits.set(2, count != 0);

					kmerFreq_t kmer_freqs   =  leaf -> fwdInterval.size() + leaf->rvcInterval.size();
					kmerFreq_t allKmer_freq =  leaf -> getKmerCount();
					double     aveKmer_freq = (double) allKmer_freq/(double) m_currentLength;
					double     kmer_ratio   = aveKmer_freq / maximumAverfreqs;

					std::string currString  =  leaf -> getSuffix(m_currentKmerSize);
					std::string debugBitStr =
							pDebugBits.to_string<char, std::string::traits_type, std::string::allocator_type>();

					size_t currLen        = totalPathInfo.size();
					size_t globalGCNumber = (lastPathInfo.at(leaf).GCNumber   );
					size_t localGCNumber  = globalGCNumber;
					size_t Homopolymer    = (lastPathInfo.at(leaf).Homopolymer);

					if (currLen >= m_currentKmerSize + 1)
					{
						size_t startLoc = currLen - m_currentKmerSize - 1;
						localGCNumber  -= (leaf -> getParentInfo(totalPathInfo,startLoc)).GCNumber;
					}

					char strand;
					if (m_extSeeds.isPosStrand)
						strand = '+';
					else
						strand = '-';

					char pruned= '-';
					if      (leaf -> removeType ==  1)
						pruned = 'O';
					else if (leaf -> removeType == -1)
						pruned = 'X';

					fieldNum = 1;
					strBuffer   << std::fixed << std::setprecision(5)
								<< ">" << fieldNum
								<< "|" << m_Debug.readID
								<< "|" << m_extSeeds.source.start
								<< "|" << m_extSeeds.source.end
								<< "|" << m_extSeeds.target.start
								<< "|" << m_extSeeds.target.end
								<< "|" << strand;

					fieldNum++;
					strBuffer   << "&" << fieldNum
								<< "|" << m_extSeeds.source.isRepeat
								<< "|" << m_extSeeds.target.isRepeat
								<< "|" << m_disBetweenSrcTarget;

					fieldNum++;
					strBuffer   << "&" << fieldNum
								<< "|" << m_case_number
								<< "|" << m_step_number;

					fieldNum++;
					strBuffer   << "&" << fieldNum
								<< "|" << m_currLeavesTotalSize
								<< "|" << m_currLeavesTotalSize - absolutePruning
								<< "|" << m_currLeavesTotalSize -(absolutePruning + relativePruning);

					fieldNum++;
					strBuffer   << "&" << fieldNum
								<< "|" << (leaf -> currLeafID)
								<< "|" << (leaf -> lastLeafID);

					fieldNum++;
					strBuffer   << "&" << fieldNum
								<< "|" << m_accumLeaves
								<< "|" << m_totalExtNum;

					fieldNum++;
					strBuffer   << "&" << fieldNum
								<< "|" << m_currentKmerSize
								<< "|" << m_currentLength;

					fieldNum++;
					strBuffer   << "&" << fieldNum
								<< "|" << kmer_freqs
								<< "|" << allKmer_freq
								<< "|" << kmer_ratio;

					fieldNum++;
					strBuffer   << "&" << fieldNum
								<< "|" << localGCNumber
								<< "|" << globalGCNumber
								<< "|" << Homopolymer;

					fieldNum++;
					strBuffer   << "&" << fieldNum
								<< "|" << 100*(lastPathInfo.at(leaf).localErrorRate )
								<< "|" << 100*(lastPathInfo.at(leaf).globalErrorRate);

					fieldNum++;
					strBuffer   << "&" << fieldNum
								<< "|" << (leaf -> totalMatch) + m_matchedSize-1
								<< "|" << (leaf -> numRedeemMatched)
								<< "|" << (leaf -> lastMatchedIdxOffset)
								<< "|" << (leaf -> lastOverlapLen)
								<< "|" << (leaf -> queryOverlapLen);

					fieldNum++;
					strBuffer   << "&" << fieldNum
								<< "|" << debugBitStr
								<< "|" << pruned;

					fieldNum++;
					strBuffer   << "&" << fieldNum
								<< "|";

				}

			if (leaf -> isKept())
			{
				// extension
					extensions = getFMIndexExtensions(leaf, pDebugBits, strBuffer);

				// Set the debug strings
					if (m_Debug.isDebug)
					{
						strBuffer << leaf -> getFullString() << "\n";
					}

				// Update the nodes after successful extension .
					if (extensions.size() > 0)
					{
						updateLeaves(extensions, leaf, lastPathInfo, currPathInfo, currLeavesNum);
						break;
					}

				// Relax the leaf with the minimum error rate.
					m_min_SA_threshold--;
					count++;
			}
			else
			{
				// Set the debug strings
					if (m_Debug.isDebug)
					{
						debugPerExtArray tempDebugData;
						tempDebugData.resize(4);
						printDebugData(tempDebugData, strBuffer);
						strBuffer << "0|0\n" << leaf -> getFullString() << "\n";
					}
				break;
			}

		}

		// Recover the threshold already relaxed.
			m_min_SA_threshold += count;

		// Compute the debug parameters
			m_currTotalExtNum  += m_currExtNum;

			if (minTotalcount >= totalcount)
				minTotalcount  = totalcount;

			if (m_Debug.isDebug)
				debugStr += strBuffer.str();

	}

	// Set the number of the next leaves
		m_nextLeavesKeptSize = m_nextLeavesTotalSize = currLeavesNum;
}

// Update the leaves after the extension is successful.
void LongReadSelfCorrectByFMExtend::updateLeaves(extArray& extensions,SAIOverlapNode3* leaf, const leafTable_t& lastPathInfo, leafTable_t& currPathInfo, size_t& currLeavesNum)
{
	const pathInfo_t& lastLeafInfo = lastPathInfo.at(leaf);
	if(extensions.size() == 1)
	{
		// Single extension, do not branch
			leaf -> extend (extensions.front().SearchLetters );

		// Set the current Leaf ID.
			currLeavesNum++;
		// Update the leaves
			pathInfo_t  pathInfo( lastLeafInfo, extensions.front() );
			leaf -> setNode(extensions.front(), currLeavesNum, leaf );
			m_nextLeaves.emplace_back(leaf);
			currPathInfo.emplace( leaf, pathInfo );

	}
	else if(extensions.size() > 1)
	{
		// Branch
		for(size_t i = 0; i < extensions.size(); ++i)
		{
			SAIOverlapNode3* pChildNode = leaf -> createChild(extensions.at(i).SearchLetters);
			//inherit accumulated kmerCount from parent
				pChildNode -> addKmerCount( leaf -> getKmerCount() );
/*
			// Remove the bug while changing the seed by FM-Extend
				if (i != 0)
					(pChildNode -> resultindex).first = -1;
//*/
			// Set the current Leaf ID.
				currLeavesNum++;
			// Update the leaves
				pathInfo_t  pathInfo( lastLeafInfo, extensions.at(i));
				pChildNode -> setNode(extensions.at(i), currLeavesNum, leaf, totalPathInfo.size());
				m_nextLeaves.emplace_back(pChildNode);
				currPathInfo.emplace( pChildNode, pathInfo );
		}

		m_accumLeaves--;
		m_accumLeaves += extensions.size();
	}
}

// Compute the error rates and keep the leaves whose error rates <= expected error rate.
bool LongReadSelfCorrectByFMExtend::PrunedByAbsErrorRate()
{
	// the matched index in m_TerminatedIntervals for m_currentLength
	// the m_currentLength is the same for all leaves
	// which is used as the central index within the m_maxIndelSize window
	size_t currMatchedIdx = m_currentLength-m_matchedSize;

	//        ---	matched size =3. matched dist = 1;
	//         --*
	//          -*-
	//           *--
	//            ---
	// *: SNP or indel errors, ---: matched size
	// Erase the leaf if no feasible matches are found within matchedSize + maxIndelSize.
	size_t indelOffset = m_matchedSize+m_maxIndelSize;

	// Compute the range of small and large indices for tolerating m_maxIndelSize
	size_t smallMatchedIdx = currMatchedIdx <= indelOffset ? 0 : currMatchedIdx - indelOffset;
	size_t largeMatchedIdx =(currMatchedIdx+indelOffset) >= (m_query.length()-m_matchedSize)?
							(m_query.length()-m_matchedSize):currMatchedIdx+indelOffset;

	// check range of last matched and find new matches for each interval
	for( auto leaf: m_nextLeaves )
	{
		bool isNewMatchFound = false;

		if (m_currentLength - leaf -> lastOverlapLen > m_matchedSize
		||  m_currentLength - leaf -> lastOverlapLen <= 1 )
		{
			size_t preMatchedIdx = leaf -> lastMatchedIdx;
			// search for a new match
			isNewMatchFound = findNewMatch(leaf, smallMatchedIdx, largeMatchedIdx);

			// lastMatchedIdxOffset records the offset between lastMatchedIdx and currMatchedIdx when first match is found
			if(isNewMatchFound)
			{
				if( currMatchedIdx + leaf -> lastMatchedIdxOffset - preMatchedIdx > m_matchedSize )
					leaf -> numRedeemMatched += (m_matchedSize-1)*m_PacBioErrorRate;

				leaf -> lastMatchedIdxOffset = (int) leaf->lastMatchedIdx - (int)currMatchedIdx;
			}
			else
			{
				// If the matched extension is stopped by SNP or indel error for the 1st time
				// increment the error number in order to distinguish two separate matches
				// and one larger consecutive matched during error rate computation
				if     ( (currMatchedIdx + leaf->lastMatchedIdxOffset - leaf->lastMatchedIdx) % m_matchedSize == 1 )
					leaf->numAscertainMismatched ++;
				else if( (currMatchedIdx + leaf->lastMatchedIdxOffset - leaf->lastMatchedIdx) > m_matchedSize  - 1 )
					leaf->numRedeemMatched += 1-m_PacBioErrorRate;
			}
		}

		else
			leaf->numRedeemMatched += 1-m_PacBioErrorRate;

		double currErrorRate = computeErrorRate(leaf);

		// speedup by skipping dissimilar reads
		// This is the 2nd filter less reliable than the 1st one
		if((!m_Debug.ref.isSpecifiedPath) && (currErrorRate > m_errorRate)) //testcw
		{
			leaf -> removeType = 1;
			m_nextLeavesKeptSize--;
		}

	}

	return true;
}

// Find the matched k-mer with the current extension sequence in the raw read.
bool LongReadSelfCorrectByFMExtend::findNewMatch(SAIOverlapNode3* currNode, size_t smallMatchedIdx, size_t largeMatchedIdx)
{
	// If there is mismatch/indel, jump to the next m_matchedSize/m_matchedDist, and 1 otherwise.
	size_t matchedIdxOffset =   currNode->lastOverlapLen < m_currentLength-m_matchedSize?
								m_matchedSize:m_currentLength - currNode->lastOverlapLen;

	// search for new matched starting from last matched or smallMatchedIdx
	size_t startMatchedIdx  = std::max(smallMatchedIdx, currNode->lastMatchedIdx+matchedIdxOffset);

	bool isNewMatchFound = false;
	BWTInterval currFwdInterval = currNode->fwdInterval;
	BWTInterval currRvcInterval = currNode->rvcInterval;

	// Binary search for new matches using Query interval tree
	std::vector<TreeInterval<size_t> > resultsFwd, resultsRvc;
	if(currFwdInterval.isValid())
		m_fwdIntervalTree.findOverlapping(currFwdInterval.lower, currFwdInterval.upper, resultsFwd);
	if(currRvcInterval.isValid())
		m_rvcIntervalTree.findOverlapping(currRvcInterval.lower, currRvcInterval.upper, resultsRvc);
	int minIdxDiff = 10000;
	size_t currMatchedIdx = m_currentLength-m_matchedSize;
	for(size_t i=0 ; i<resultsFwd.size() || i<resultsRvc.size() ; i++)
	{
		if(
				currFwdInterval.isValid()
		&&      i < resultsFwd.size()
		&&      resultsFwd.at(i).value >= startMatchedIdx
		&&      resultsFwd.at(i).value <= largeMatchedIdx
			)
		{
			// update currNode members
			if(std::abs((int)resultsFwd.at(i).value - (int)currMatchedIdx) < minIdxDiff)
			{
				currNode->lastMatchedIdx = resultsFwd.at(i).value;

				// query overlap may shift due to indels
					currNode->queryOverlapLen = resultsFwd.at(i).value+m_matchedSize;
					minIdxDiff = std::abs((int)resultsFwd.at(i).value - (int)currMatchedIdx);
			}
			// lastOverlapLen records the overlap length of last hit
				currNode->lastOverlapLen = m_currentLength;
			// currOverlapLen is always identical to m_currentLength
				currNode->currOverlapLen = m_currentLength;
			isNewMatchFound = true;
		}
		else if(
					currRvcInterval.isValid()
			&&      i < resultsRvc.size()
			&&      resultsRvc.at(i).value >= startMatchedIdx
			&&      resultsRvc.at(i).value <= largeMatchedIdx
				)
		{

			// update currNode members
			if(std::abs( (int)currMatchedIdx - (int)resultsRvc.at(i).value ) < minIdxDiff)
			{
				currNode->lastMatchedIdx = resultsRvc.at(i).value;

				// query overlap may shift due to indels
					currNode->queryOverlapLen = resultsRvc.at(i).value+m_matchedSize;
					minIdxDiff = std::abs((int)resultsRvc.at(i).value - (int)currMatchedIdx);
			}
			// lastOverlapLen records the overlap length of last hit
				currNode->lastOverlapLen = m_currentLength;
			// currOverlapLen is always identical to m_currentLength
				currNode->currOverlapLen = m_currentLength;
			isNewMatchFound = true;
		}
	}

	if(isNewMatchFound)
		currNode->totalMatch++;

	return isNewMatchFound;
}

// Compute the error rate in the leaf.
double LongReadSelfCorrectByFMExtend::computeErrorRate(SAIOverlapNode3* currNode)
{
	// Compute accuracy via matched length in both query and subject

	double matchedLen = (double)currNode->totalMatch + m_matchedSize-1;

	// SNP and indel over-estimate the unmatched lengths across error, ---*---
	// Restore the unmatched region via numAscertainMismatched, which is still over-estimated

	matchedLen += currNode->numRedeemMatched;

	double totalLen = (double)currNode->currOverlapLen ;

	double unmatchedLen = totalLen - matchedLen;

	double  globalErrorRate = unmatchedLen/totalLen;
	double  localErrorRate  = globalErrorRate;

	// Set the offset mode
	size_t totalsize = totalPathInfo.size();
	if( totalsize >= m_localSimilarlykmerSize + m_pathOffset)
	{
		size_t       startLoc       = totalsize - m_localSimilarlykmerSize;
		size_t       startLength    = totalLen  - m_localSimilarlykmerSize;
		const double startErrorRate =(currNode  -> getParentInfo(totalPathInfo, startLoc)).globalErrorRate;

		localErrorRate = ( globalErrorRate * totalLen - startErrorRate * startLength )/m_localSimilarlykmerSize;
	}

	(currNode -> getParentInfo(totalPathInfo)).setErrorRate( localErrorRate, globalErrorRate);

	return localErrorRate;
}

// Attempt to the extension in the current leaf and get the extension information.
extArray LongReadSelfCorrectByFMExtend::getFMIndexExtensions(SAIOverlapNode3* leaf, const debugBits_t& debugBits, std::stringstream& strBuffer)
{
// Set the array data which is reserved.
	extArray output;
		output.reserve(4);
	extArray totalExt;
		totalExt.reserve(4);
	debugPerExtArray debugData;
		debugData.reserve(4);
	std::vector<bool> extMatch;
		extMatch.reserve(4);

	size_t IntervalSizeCutoff = m_min_SA_threshold;
						//min freq at fwd and rvc bwt, >=3 is equal to >=2 kmer freq
	int maxfreqsofleaf = 0;
	totalcount = 0;
	m_currExtNum = 0;
	dominantBase maxFreqOfLeaves(leaf -> getSuffix(1));

	for(int i = 1; i < BWT_ALPHABET::size; ++i) //i = A,C,G,T
	{
		//update forward Interval using extension b
			char b = BWT_ALPHABET::getChar(i);
			BWTInterval fwdInterval = leaf->fwdInterval;
			if(fwdInterval.isValid())
				BWTAlgorithms::updateInterval(fwdInterval,b,m_pRBWT);

		//update reverse complement Interval using extension rcb
			char rcb=BWT_ALPHABET::getChar(5-i);
			BWTInterval rvcInterval = leaf->rvcInterval;
			if(rvcInterval.isValid())
				BWTAlgorithms::updateInterval(rvcInterval,rcb,m_pBWT);

		FMidx_t currExt       = FMidx_t(b,fwdInterval,rvcInterval);
		bool    isMatchedLast = ismatchedbykmer(fwdInterval,rvcInterval);
		kmerFreq_t currFreqs  = currExt.getKmerFrequency();

		totalcount += currFreqs;

		if (currFreqs != 0)
			m_currExtNum++;

		maxFreqOfLeaves.setFreq(currExt);
		totalExt.push_back(currExt);
		extMatch.push_back(isMatchedLast);
		debugData.emplace_back(currExt, isMatchedLast);

	}// end of ACGT

	m_maxfreqs = std::max(m_maxfreqs,totalcount);
	maxfreqsofleaf = maxFreqOfLeaves.getMaxFreq();

	for(int i = 1; i < BWT_ALPHABET::size; ++i)
	{
		// Get the information
		size_t kmerFreq          = totalExt.at(i-1).getKmerFrequency();
		BWTInterval fwdInterval  = totalExt.at(i-1).getFwdInterval();
		BWTInterval rvcInterval  = totalExt.at(i-1).getRvcInterval();
		const bool isMatchedLast = extMatch.at(i-1);

		// Compute the k-mer ratio arguments
		const double kmerRatioNotPass = 2;
		double kmerRatioCutoff = 0;
		double kmerRatio = (double) kmerFreq/(double)maxfreqsofleaf;

		char b = BWT_ALPHABET::getChar(i);

		bool isHomopolymer = ((leaf -> getParentInfo(totalPathInfo)).Homopolymer >= 3);

		bool isDominant     = maxFreqOfLeaves.isDominant();
		bool isFreqPass     = kmerFreq   >= IntervalSizeCutoff;
		bool isLowCoverage  = totalcount >= IntervalSizeCutoff+2;
		bool isRepeat       = maxfreqsofleaf > 100;
		bool isHighlyRepeat = maxfreqsofleaf > 150;
		bool isLowlyRepeat  = maxfreqsofleaf >  50;

		// matched case (repeat case)
			if  ( isMatchedLast &&  isHighlyRepeat )
			{
				kmerRatioCutoff = 0.125;
				debugData.at(i-1).setErrorData(kmerRatio, "MHR");   // Match High Repeat
			}
		else if ( isMatchedLast &&  isLowlyRepeat  )
			{
				kmerRatioCutoff = 0.2  ;
				debugData.at(i-1).setErrorData(kmerRatio, "MLR");   // Match Low Repeat
			}
		// unmatched case (unique case)
		else if ( isFreqPass )
			{
				kmerRatioCutoff = 0.25 ;
				debugData.at(i-1).setErrorData(kmerRatio, "NCU");   // Normal Coverage Unique
			}
		else if ( isLowCoverage )
			{
				kmerRatioCutoff = 0.6  ;
				debugData.at(i-1).setErrorData(kmerRatio, "NLC");   // Normal Low Coverage
			}
		else
			{
				kmerRatioCutoff = kmerRatioNotPass;
				debugData.at(i-1).setErrorData(kmerRatio, "TLC");   // Too Low Coverage
			}

		// Homopolymer case
			if  ( isHomopolymer   &&  isRepeat )
			{
				if (kmerRatioCutoff <= 0.3)
				{
					kmerRatioCutoff = 0.3;
					debugData.at(i-1).setErrorData(kmerRatio, "HRR"); // Homopolymer Repeat Region
				}
			}
		else if ( isHomopolymer )
			{
				if (kmerRatioCutoff <= 0.6)
				{
					kmerRatioCutoff = 0.6;
					debugData.at(i-1).setErrorData(kmerRatio, "HUR"); // Homopolymer Unique Region
				}
			}

		if( kmerRatio >=   kmerRatioCutoff )
		{
			// extend to b
				output.emplace_back(b,fwdInterval,rvcInterval);
				debugData.at(i-1).extSuccessSet();
		}
	}

	// The debug mode to force to extend into the correct path.
		if(m_Debug.ref.isSpecifiedPath)
		{
			if (!(m_step_number == 1 || m_Debug.ref.isMatchBase(leaf -> getSuffix(1), m_step_number-2)))
			{
				output.clear();
			}
			else if(!output.empty() || debugBits[1] || debugBits[2])
			{
				// output is empty, and the finally thresholds reduced
				// output is not empty
				output.clear();
				for(int i = 1; i < BWT_ALPHABET::size; ++i)
				{
					FMidx_t&           currExt     = totalExt.at(i-1);
					const BWTInterval& fwdInterval = currExt.getFwdInterval();
					const BWTInterval& rvcInterval = currExt.getRvcInterval();
					char b = BWT_ALPHABET::getChar(i);

					if  (m_step_number == 1)
						{
							if (currExt.isVaildIntervals())
								output.emplace_back(b,fwdInterval,rvcInterval);
						}
					else if  (currExt.isVaildIntervals() || m_Debug.ref.isMatchBase(b, m_step_number-1))
						{
							output.emplace_back(b,fwdInterval,rvcInterval);
						}
				}
			}
		}

	// Set the debug strings
		if(m_Debug.isDebug)
		{
			printDebugData(debugData, strBuffer);
			strBuffer   << maxfreqsofleaf << "|"
						<< totalcount     << "\n";
		}

	return output;
}

// Determine if the current sequence is matched 5-mer raw read.
bool LongReadSelfCorrectByFMExtend::ismatchedbykmer(BWTInterval currFwdInterval,BWTInterval currRvcInterval)
{
	bool match = false;

	// Binary search for new matches using Query interval tree
	std::vector<TreeInterval<size_t> > resultsFwd, resultsRvc;
	if(currFwdInterval.isValid())
		m_fwdIntervalTree2.findOverlapping(currFwdInterval.lower, currFwdInterval.upper, resultsFwd);
	if(currRvcInterval.isValid())
		m_rvcIntervalTree2.findOverlapping(currRvcInterval.lower, currRvcInterval.upper, resultsRvc);
	size_t startMatchedIdx = std::max((int)m_currentLength - (int)m_maxIndelSize,0);
	size_t largeMatchedIdx = m_currentLength + m_maxIndelSize;

	for(size_t i=0 ; i<resultsFwd.size() || i<resultsRvc.size() ; i++)
	{
		if(
			currFwdInterval.isValid()
		&&  i < resultsFwd.size()
		&&  resultsFwd.at(i).value >= startMatchedIdx
		&&  resultsFwd.at(i).value <= largeMatchedIdx
			)
		{
			match = true;
			break;
		}
		else if(
				currRvcInterval.isValid()
			&&  i < resultsRvc.size()
			&& resultsRvc.at(i).value >= startMatchedIdx
			&& resultsRvc.at(i).value <= largeMatchedIdx
			)
		{
			match = true;
			break;
		}
	}

	return match;
}

// Check for leaves whose extension has terminated.
// If the leaf has terminated, the walked string and coverage is pushed to the result vector
bool LongReadSelfCorrectByFMExtend::isTerminated(SAIntervalNodeResultVector& results)
{
	bool found = false;

	for(auto leaf : m_currLeaves)
	{
		if (leaf -> isRemoved())
			continue;

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
				STresult.errorRate   = (leaf->getParentInfo(totalPathInfo)).globalErrorRate;
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

void LongReadSelfCorrectByFMExtend::printDebugData(debugPerExtArray& currDebugData, std::stringstream& strBuffer)
{
// Print the error type.
	for(int i = 1; i < BWT_ALPHABET::size; ++i)
	{
		strBuffer << currDebugData.at(i-1).errorType;
		if (i == BWT_ALPHABET::size-1)
		{
			fieldNum++;
			strBuffer   << "&" <<  fieldNum ;
		}
		strBuffer << "|" ;
	}
// Print the k-mer frequency.
	for(int i = 1; i < BWT_ALPHABET::size; ++i)
	{
		strBuffer << currDebugData.at(i-1).kmerFreq;
		if (i == BWT_ALPHABET::size-1)
		{
			fieldNum++;
			strBuffer   << "&" <<  fieldNum ;
		}
		strBuffer << "|" ;
	}
// Print the k-mer ratio.
	for(int i = 1; i < BWT_ALPHABET::size; ++i)
	{
		strBuffer << currDebugData.at(i-1).kmerRatio;
		if (i == BWT_ALPHABET::size-1)
		{
			fieldNum++;
			strBuffer   << "&" <<  fieldNum ;
		}
		strBuffer << "|" ;
	}
// Print "is matched".
	for(int i = 1; i < BWT_ALPHABET::size; ++i)
	{
		strBuffer << currDebugData.at(i-1).isMatched;
		if (i == BWT_ALPHABET::size-1)
		{
			fieldNum++;
			strBuffer   << "&" <<  fieldNum ;
		}
		strBuffer << "|" ;
	}

}

void LongReadSelfCorrectByFMExtend::printErrorRate(LeavesArray_t& currLeaves)
{
	std::cout << std::fixed << std::setprecision(5);

	for(auto& leaf : currLeaves)
	{
		if ( leaf -> isRemoved() )
			continue;

		std::cout << 100*((leaf->getParentInfo(totalPathInfo)).localErrorRate) << "\t";
	}

	std::cout << std::endl;
}