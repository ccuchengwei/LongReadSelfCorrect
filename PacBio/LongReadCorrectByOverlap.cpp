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
					m_Debug(debug),
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
		freqsOfKmerSize = new double[100 + 1]{0};
		for(int i = m_minOverlap ; i <= 100 ; i++)
			freqsOfKmerSize[i] = pow(1 - m_PacBioErrorRate, i) * m_PBcoverage;

	// PacBio reads are longer than real length due to insertions
		m_maxLength =          (1.2*(m_disBetweenSrcTarget+10))+2*m_initkmersize     ;
		m_minLength = std::max((0.8*(m_disBetweenSrcTarget-20))+2*m_initkmersize,0.0);

	//initialRootNode
		initialRootNode(m_extSeeds.source.seq);

	// push new node into roots and leaves vector
		m_RootNodes.push_back(m_pRootNode);
		m_leaves.emplace_back(m_pRootNode);

		totalPathInfo.reserve(m_maxLength+1);
		pathInfo_t  pathInfo(m_extSeeds.source.seq);
		totalPathInfo.emplace_back
			(
				leafTable_t ({
								{m_pRootNode, pathInfo}
							})
			);
//*
	if(m_Debug.isDebug)
	{
		(*m_Debug.debug_finalSeq)
					<< m_Debug.caseNum << "\tBE: "
					<< m_extSeeds.source.seq << " "
					<< m_extSeeds.target.seq << " "
					<< m_pRootNode->fwdInterval.size() + m_pRootNode->rvcInterval.size() << "|"
					<< disBetweenSrcTarget <<"\n";
	}
//*/

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
		m_pRootNode -> setRoot(beginningkmer, m_initkmersize, m_seedSize, m_pBWT, m_pRBWT);

		m_maxfreqs = m_pRootNode->fwdInterval.size() + m_pRootNode->rvcInterval.size();
		m_currentLength = m_currentKmerSize = m_initkmersize;
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
			SONode3PtrList newLeaves;
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
	size_t realCandiNum = 1;
	size_t mismatchSNum = 1;

	for (size_t i = 0 ; i < results.size() ;i++)
	{
		const std::string& candidateSeq = results[i].thread;
//*
		if(m_Debug.isDebug)
		{
			(*m_Debug.debug_finalSeq)
							<< "Final Seqs: "  << results[i].thread
							<< "\tError Rate:" << results[i].errorRate << std::endl;
		}
//*/
		if( results[i].errorRate < minErrorRate )
		{
			realCandiNum = 1;
			mismatchSNum = 1;
			minErrorRate = results[i].errorRate;
			FMWResult.mergedSeq = candidateSeq;
			minTotalcount = results[i].SAIntervalSize;
		}
		else if (results[i].errorRate == minErrorRate )
		{
			realCandiNum++;
			if (m_Debug.isDebug)
			{
				if (candidateSeq.length() == FMWResult.mergedSeq.length())
				{
					mismatchSNum++;
				}
			}
		}
	}
/*
	if (m_Debug.isDebug)
	{
		(*m_Debug.debug_finalSeq)
						<< m_Debug.readID        << "\t"
						<< m_Debug.caseNum       << "\t"
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
void LongReadSelfCorrectByOverlap::extendLeaves(SONode3PtrList& newLeaves)
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

	m_totalExtNum += m_currTotalExtNum;

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
size_t LongReadSelfCorrectByOverlap::SelectFreqsOfrange(const size_t LowerBound, const size_t UpperBound, SONode3PtrList& newLeaves)
{
	extArray maxKmerArray;
	int tempmaxfmfreqs = 0;

	for(auto& leaf : newLeaves)
	{
		std::string    maxKmer = leaf -> getSuffix(UpperBound);

		std::string startkmer  = maxKmer.substr(UpperBound - LowerBound); //  string of lower bound kmer size

		BWTInterval Fwdinterval = BWTAlgorithms::findInterval(m_pBWT, startkmer);
		BWTInterval Rvcinterval = BWTAlgorithms::findInterval(m_pRBWT, reverseComplement(reverse(startkmer)));

		maxKmerArray.emplace_back(maxKmer,Fwdinterval,Rvcinterval);
		FMidx_t& currKmer = maxKmerArray.back();

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
bool LongReadSelfCorrectByOverlap::isInsufficientFreqs(SONode3PtrList& newLeaves)
{
	size_t highfreqscount = 0;
	for(auto& leaf : newLeaves)
	{
		size_t highfreqThreshold = m_PBcoverage > 60 ? (size_t)(m_PBcoverage/60)*3 : 3;

		if( (leaf -> getLastKmerCount()) > highfreqThreshold)
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
void LongReadSelfCorrectByOverlap::refineSAInterval(SONode3PtrList& leaves, const size_t newKmerSize)
{
	for(auto& leaf : leaves)
	{
		// reset the SA intervals using newKmerSize
			std::string reducedKmer = leaf -> getSuffix(newKmerSize);

			leaf -> fwdInterval = BWTAlgorithms::findInterval(m_pRBWT,           reverse(reducedKmer));
			leaf -> rvcInterval = BWTAlgorithms::findInterval(m_pBWT , reverseComplement(reducedKmer));
	}

	m_currentKmerSize = newKmerSize;
}

// Keep the leaves whose highly-correct rates are relative to the others.
// And attempt to extend those leaves.
void LongReadSelfCorrectByOverlap::attempToExtend(SONode3PtrList& newLeaves, const bool isSuccessToReduce)
{

	const   leafTable_t& lastPathInfo = totalPathInfo.back();
			leafTable_t  currPathInfo;     // the max branches = 4 per leaf.
				currPathInfo.reserve(4 * lastPathInfo.size());

	m_maxfreqs = 0;

	std::vector<size_t> frequencies;


	// Compute the min error rate
		double minimumErrorRate = 1.0;
		for(auto& leaf : m_leaves)
		{
			if( lastPathInfo.at(leaf).localErrorRate < minimumErrorRate)
				minimumErrorRate = lastPathInfo.at(leaf).localErrorRate;
		}

	// Compute the errorRateDiff to trim leaves whose error rates relative to the others is high.
		SONode3PtrList::iterator iter = m_leaves.begin();
		while(iter != m_leaves.end())
		{
			if (m_Debug.ref.isSpecifiedPath)
				break;

			SAIOverlapNode3* leaf = (*iter);
			double errorRateDiff  = lastPathInfo.at(leaf).localErrorRate - minimumErrorRate;
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
	m_currTotalExtNum  = 0;

	double maximumAverfreqs = 0.0;
	if (m_Debug.isDebug && isSuccessToReduce)
	{
		for(auto& leaf : m_leaves)
		{
			kmerFreq_t allKmer_freq = leaf -> getKmerCount();
			double     aveKmer_freq = (double) allKmer_freq/(double) m_currentLength;

			if( maximumAverfreqs < aveKmer_freq )
				maximumAverfreqs = aveKmer_freq ;
		}

	}

	iter = m_leaves.begin();
	while(iter != m_leaves.end())
	{
		extArray extensions;
		int count = 0;
		SAIOverlapNode3* leaf = (*iter);
		bool isReduced = isSuccessToReduce;
		while(count < 2)
		{
			if	( count == 1
			&& !(lastPathInfo.at(leaf).localErrorRate == minimumErrorRate && m_leaves.size() > 1))
				break;

			if (m_Debug.isDebug && isReduced)
			{
				kmerFreq_t kmer_freqs   =  leaf -> fwdInterval.size() + leaf->rvcInterval.size();
				kmerFreq_t allKmer_freq =  leaf -> getKmerCount();
				double     aveKmer_freq = (double) allKmer_freq/(double) m_currentLength;
				double     kmer_ratio   = aveKmer_freq / maximumAverfreqs;

				std::string currString  =  leaf -> getSuffix(m_currentKmerSize);

				size_t currLen     = totalPathInfo.size();
				size_t GCNumber    = (lastPathInfo.at(leaf).GCNumber   );
				size_t Homopolymer = (lastPathInfo.at(leaf).Homopolymer);

				if (currLen >= m_localSimilarlykmerSize)
				{
					size_t startLoc = currLen - m_localSimilarlykmerSize;
					GCNumber -= (leaf -> getParentInfo(totalPathInfo,startLoc)).GCNumber;
				}

				char strand;
				if(m_extSeeds.isPosStrand)
					strand = '+';
				else
					strand = '-';

				*(m_Debug.debug_file)   << std::fixed << std::setprecision(5)
										<< ">" << m_Debug.readID
										<< "|" << m_extSeeds.source.start
										<< "|" << m_extSeeds.source.end
										<< "|" << m_extSeeds.target.start
										<< "|" << m_extSeeds.target.end
										<< "|" << strand

										<< "&" << m_extSeeds.source.isRepeat
										<< "|" << m_extSeeds.target.isRepeat
										<< "|" << m_disBetweenSrcTarget

										<< "&" << m_Debug.caseNum
										<< "|" << m_step_number

										<< "&" << m_leaves.size()
										<< "|" << currLeavesNum
										<< "|" << (leaf -> lastLeafID)
										<< "&" << m_accumLeaves
										<< "|" << m_totalExtNum

										<< "&" << m_currentKmerSize
										<< "|" << m_currentLength

										<< "&" << GCNumber
										<< "|" << Homopolymer

										<< "&" << kmer_freqs
										<< "|" << allKmer_freq
										<< "|" << kmer_ratio

										<< "&" << 100*(lastPathInfo.at(leaf).localErrorRate )
										<< "|" << 100*(lastPathInfo.at(leaf).globalErrorRate)

										<< "&" << (leaf -> totalSeeds) + m_seedSize-1
										<< "|" << (leaf -> numRedeemSeed)
										<< "|" << (leaf -> lastSeedIdxOffset)
										<< "|" << (leaf -> lastOverlapLen)
										<< "|" << (leaf -> queryOverlapLen)
										<< "&";
			}

			extensions = getFMIndexExtensions(leaf, isReduced);

			if (m_Debug.isDebug && isReduced)
			{
				*(m_Debug.debug_file) << leaf -> getFullString() << "\n";
			}

			if(extensions.size() > 0)
			{
				updateLeaves(newLeaves, extensions, leaf, lastPathInfo, currPathInfo, currLeavesNum);
				break;
			}
			isReduced = false;
			m_min_SA_threshold--;
			count++;
		}
		m_min_SA_threshold += count;
		m_currTotalExtNum += m_currExtNum;

		if (minTotalcount >= totalcount)
		{
			minTotalcount = totalcount;
		}

		++iter;
		++currLeavesNum;
	}
	if (!newLeaves.empty())
	{
		totalPathInfo.emplace_back(currPathInfo);
	}
}

// Update the leaves after the extension is successful.
void LongReadSelfCorrectByOverlap::updateLeaves(SONode3PtrList& newLeaves,extArray& extensions,SAIOverlapNode3* leaf, const leafTable_t& lastPathInfo, leafTable_t& currPathInfo, size_t currLeavesNum)
{
	if(extensions.size() == 1)
	{
		// Single extension, do not branch
			leaf -> extend (extensions.front().SearchLetters);
			leaf -> setNode(extensions.front(), currLeavesNum, leaf);
			newLeaves.emplace_back(leaf);

			pathInfo_t  pathInfo( lastPathInfo.at(leaf), extensions.front());
			currPathInfo.emplace( leaf, pathInfo);
	}
	else if(extensions.size() > 1)
	{
		// Branch
		const pathInfo_t& lastLeafInfo = lastPathInfo.at(leaf);
		for(size_t i = 0; i < extensions.size(); ++i)
		{
			SAIOverlapNode3* pChildNode = leaf -> createChild(extensions[i].SearchLetters);
			//inherit accumulated kmerCount from parent
				pChildNode -> addKmerCount( leaf -> getKmerCount() );
			pChildNode -> setNode(extensions[i], currLeavesNum, leaf, m_step_number);
/*
			// Remove the bug while changing the seed by FM-Extend
				if (i != 0)
					(pChildNode -> resultindex).first = -1;
//*/
			newLeaves.emplace_back(pChildNode);

			pathInfo_t  pathInfo( lastLeafInfo, extensions[i]);
			currPathInfo.emplace( pChildNode, pathInfo);
		}

		m_accumLeaves--;
		m_accumLeaves += extensions.size();
	}
}

// Compute the error rates and keep the leaves whose error rates <= expected error rate.
bool LongReadSelfCorrectByOverlap::PrunedBySeedSupport(SONode3PtrList& newLeaves)
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
	SONode3PtrList::iterator iter = newLeaves.begin();
	while(iter != newLeaves.end())
	{
		bool isNewSeedFound = false;
		SAIOverlapNode3* leaf = (*iter);

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

	double  globalErrorRate = unmatchedLen/totalLen;
	double  localErrorRate  = globalErrorRate;

	if( totalPathInfo.size() >= m_localSimilarlykmerSize )
	{
		size_t       totalsize      = totalPathInfo.size();
		size_t       startLoc       = totalsize - m_localSimilarlykmerSize;
		size_t       startLength    = totalLen  - m_localSimilarlykmerSize;
		const double startErrorRate =(currNode -> getParentInfo(totalPathInfo, startLoc)).globalErrorRate;

		localErrorRate = ( globalErrorRate * totalLen - startErrorRate * startLength )/m_localSimilarlykmerSize;
	}

	(currNode -> getParentInfo(totalPathInfo)).setErrorRate( localErrorRate, globalErrorRate);

	return localErrorRate;
}

// Attempt to the extension in the current leaf and get the extension information.
extArray LongReadSelfCorrectByOverlap::getFMIndexExtensions(SAIOverlapNode3* leaf,const bool printDebugInfo)
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
/*
	if(m_Debug.isDebug)
		std::cout   << leaf->getFullString() <<" || Local Error Rate: "
					<< (leaf->getParentInfo(totalPathInfo)).localErrorRate  <<"\n"
					<< leaf->getSuffix(m_currentKmerSize) << " || k-mer freqs: "
					<< leaf->getLastKmerCount()           << std::endl;
*/
	for(int i = 1; i < BWT_ALPHABET::size; ++i) //i=A,C,G,T
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

		// matched case
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
		// unmatched case
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

	if(m_Debug.ref.isSpecifiedPath)
	{
		if (!(m_step_number == 1 || m_Debug.ref.isMatchBase(leaf -> getSuffix(1), m_step_number-2)))
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

	if(m_Debug.isDebug && printDebugInfo)
	{
		printDebugData(debugData);
		*(m_Debug.debug_file)

			<< maxfreqsofleaf << "|"
			<< totalcount     << "\n";
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

	for(auto leaf : m_leaves)
	{
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

void LongReadSelfCorrectByOverlap::printDebugData(debugPerExtArray& currDebugData)
{
// Print the error type.
	for(int i = 1; i < BWT_ALPHABET::size; ++i)
	{
		*(m_Debug.debug_file) << currDebugData.at(i-1).errorType;
		if (i != BWT_ALPHABET::size-1)
				*(m_Debug.debug_file) << "|";
			else
				*(m_Debug.debug_file) << "&";
	}
// Print the k-mer frequency.
	for(int i = 1; i < BWT_ALPHABET::size; ++i)
	{
		*(m_Debug.debug_file) << currDebugData.at(i-1).kmerFreq;
		if (i != BWT_ALPHABET::size-1)
				*(m_Debug.debug_file) << "|";
			else
				*(m_Debug.debug_file) << "&";
	}
// Print the k-mer ratio.
	for(int i = 1; i < BWT_ALPHABET::size; ++i)
	{
		*(m_Debug.debug_file) << currDebugData.at(i-1).kmerRatio;
		if (i != BWT_ALPHABET::size-1)
				*(m_Debug.debug_file) << "|";
			else
				*(m_Debug.debug_file) << "&";
	}
// Print "is matched".
	for(int i = 1; i < BWT_ALPHABET::size; ++i)
	{
		*(m_Debug.debug_file) << currDebugData.at(i-1).isMatched;
		if (i != BWT_ALPHABET::size-1)
				*(m_Debug.debug_file) << "|";
			else
				*(m_Debug.debug_file) << "&";
	}

}

void LongReadSelfCorrectByOverlap::printErrorRate(SONode3PtrList& currLeaves)
{
	std::cout << std::fixed << std::setprecision(5);

	for(auto& leaf : currLeaves)
	{
		std::cout << 100*((leaf->getParentInfo(totalPathInfo)).localErrorRate) << "\t";
	}

	std::cout << std::endl;
}