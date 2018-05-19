///----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// LongReadExtend - Iteratively construct a
// string representing a walk through an assembly graph
// matching a query sequence.
//
// The assembly graph is abstractly represented as
// an FM-index.
//
#include "LongReadExtend.h"
#include "BWTAlgorithms.h"


//
// Class: LongReadExtend
LongReadExtend::LongReadExtend(
		const std::string& _source,
		const std::string& _target,
		int _minOverlap,
		int _maxOverlap,
		int _maxLength,
		int _maxLeaves,
		BWTIndexSet _indices,
		int _SA_threshold)
:	source(_source),
	target(_target),
	minOverlap(_minOverlap),
	maxOverlap(_maxOverlap),
	maxLength(_maxLength),
	maxLeaves(_maxLeaves),
	indices(_indices), 
	min_SA_threshold(_SA_threshold),
//	maxKmerCoverage(0),
//	maxUsedLeaves(0), 
	isBubbleCollapsed(false)
{
	// Create the root node containing the seed string
	pRootNode = new SAIPBNode(&source, nullptr);
	pRootNode->computeInitial(source);   //store initial str of root
	leaves.push_back(pRootNode);

	currentLength = source.length();
	currentKmerSize = minOverlap;

	//src kmer is a suffix of first read
	//initialize the src SA intervals with kmer length=minOverlap
	std::string srcKmer = source.substr(currentLength - minOverlap);
	pRootNode->biInterval = BWTAlgorithms::findBiInterval(indices, srcKmer);

	//trg kmer is a prefix of second read
	//initialize the trg SA intervals with kmer length=minOverlap
	std::string trgKmer = target.substr(0, minOverlap);
	destInterval = BWTAlgorithms::findBiInterval(indices, trgKmer);
}
//
LongReadExtend::~LongReadExtend()
{
	// Recursively destroy the tree
	delete pRootNode;
}

//On success return the length of merged string
int LongReadExtend::bindTwoSeeds(std::string &mergedseq)
{
	std::vector<std::string> results;
	
//	if( isTwoReadsOverlap(mergedseq))
//		return 1;

	//BFS search from 1st to 2nd seed via FM-index walk
	while(!leaves.empty() && (int)leaves.size() <= maxLeaves && currentLength <= maxLength)
	{
		// ACGT-extend the leaf nodes via updating existing SA interval
		extendLeaves();
	//	if(leaves.size() > maxUsedLeaves) maxUsedLeaves=leaves.size();
		//see if terminating string is reached
		if(isTerminated(results)) break;		
	}
		
	//find the path with maximum kmer coverage
	if(results.size() > 0)
	{
		/*
		//if multiple paths are bubbles collapsing all together at terminal
		if(results.size() == leaves.size()) isBubbleCollapsed=true;
		std::string tmpseq;
		for (size_t i = 0 ; i < results.size() ;i++){
			//bug fix: secondread may be shorter than minOverlap
			if(secondread.length()>minOverlap)
				tmpseq=results[i].thread+secondread.substr(minOverlap);
			else
				tmpseq=results[i].thread;
			
			size_t cov = calculateKmerCoverage (tmpseq, minOverlap, indices.pBWT);
			// size_t cov=results[i].SAICoverage;
			if (cov > maxKmerCoverage)
			{
				mergedseq=tmpseq;
				maxKmerCoverage=cov;
			}			
		}
		*/
		return 1;
	}
	
	//Did not reach the terminal kmer
	if(leaves.empty())
		return -1;	//high error
	else if(currentLength > maxLength)
		return -2;	//exceed search depth
	else if((int)leaves.size() > maxLeaves)
		return -3;	//too much repeats
	else
		return -4;
}

void LongReadExtend::extendLeaves()
{
	SAINodePtrList newLeaves;
	
	//attempt to extend one base for each leave
	attempToExtend(newLeaves);
	
	//shrink the SAIntervals in case overlap is larger than read length, which lead to empty newLeaves
	if(newLeaves.empty())
	{
		refineSAInterval(minOverlap);
		attempToExtend(newLeaves);
	}	
	
	//extension succeed
	if(!newLeaves.empty())
	{
		currentKmerSize++;
		currentLength++;  
	}

	leaves.clear();
	leaves = newLeaves;

	if(!leaves.empty() && currentKmerSize >= maxOverlap)
		refineSAInterval(minOverlap);

}

// Extend each leaf node
void LongReadExtend::attempToExtend(SAINodePtrList &newLeaves)
{
	for(Helper<SAINode, SAIPBNode> iter : leaves)
	{
		std::vector<std::pair<std::string, BiBWTInterval> > extensions;
		extensions = getFMIndexExtensions(iter);

		// Either extend the current node or branch it
		// If no extension, do nothing and this node
		// is no longer considered a leaf
		if(extensions.size() == 1)
		{
			// Single extension, do not branch
			const auto& probe = extensions.back();
			iter->extend(probe.first);
			iter->biInterval = probe.second;
			newLeaves.push_back(iter);
		}
		else
		{
			// Branch
			for(const auto& probe : extensions)
			{
				SAIPBNode* pChildNode = iter->createChild<SAIPBNode>(probe.first);
				pChildNode->biInterval = probe.second;
				newLeaves.push_back(pChildNode);
			}
		}
	}	
}

//update SA intervals of each leaf, which corresponds to one-base extension
std::vector<std::pair<std::string, BiBWTInterval> > LongReadExtend::getFMIndexExtensions(SAIPBNode* pNode)
{
	std::vector<std::pair<std::string, BiBWTInterval> > out;

	for(int i = 0; i < DNA_ALPHABET::size; i++) //i=A,C,G,T
	{
		char b = DNA_ALPHABET::getBase(i);
		BiBWTInterval probe = pNode->biInterval;
		BWTAlgorithms::updateBiInterval(probe, b, indices);
		
		int bcount = probe.getFreq();
		if(bcount >= min_SA_threshold)
			out.push_back(std::make_pair(std::string(1, b), probe));
	}// end of ACGT

	return out;
}

// Refine SA intervals of each leave with a new kmer
void LongReadExtend::refineSAInterval(int newKmerSize)
{
	assert(currentLength >= newKmerSize);

	for(Helper<SAINode, SAIPBNode> iter : leaves)
	{
		// reset the SA intervals using original minOverlap
		std::string newKmer = iter->getSuffix(newKmerSize);
		iter->biInterval = BWTAlgorithms::findBiInterval(indices, newKmer);
	}

	currentKmerSize = newKmerSize;
}

// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector
bool LongReadExtend::isTerminated(std::vector<std::string>& results)
{
	bool found = false;

	for(Helper<SAINode, SAIPBNode> iter : leaves)
	{
		if(iter->biInterval.isOverlapping(destInterval))
		{
			found =  found || true;
			results.push_back(iter->getFullString());
		}
	}

	return found;
}
/*
bool LongReadExtend::isTwoReadsOverlap(std::string & mergedseq)
{
	//case 1: 1st read sense overlap to 2nd read at exact minOverlap bases
	if(BWTInterval::equal(pRootNode->fwdInterval, fwdTerminatedInterval))
	{
		mergedseq= (*pQuery)+secondread.substr(minOverlap);
		return true;
	}

	//case 2: 1st read sense overlap 2nd read
	std::string secondLeftKmer=secondread.substr(0,minOverlap);
	//assume overlap can't exceed 100 bp
	size_t pos=pQuery->find(secondLeftKmer, pQuery->length()>=200?pQuery->length()-200:0);
	if(pos!=std::string::npos)	
	{
		//make sure entire suffix of 1st read after pos matches the prefix of 2nd read
		if( pQuery->substr(pos) == secondread.substr(0, pQuery->length()-pos) )
		{
			mergedseq=pQuery->substr(0,pos)+secondread;
			return true;
		}
	}

	//case 3: 1st read antisense overlap with 2nd read, or 1st read is substr of 2nd read
	//This is rare case and we don't do this in kmerMode during island joint
	if(kmerMode) return false;
	std::string firstLeftKmer=pQuery->substr(0,minOverlap);
	pos=secondread.find(firstLeftKmer);
	//assume antisense overlap can't exceed 50bp due to rare cases
	if(pos!=std::string::npos && pos <=50)
	{
		//make sure entire suffix of 2nd read after pos matches the prefix of 1st read
		if( secondread.substr(pos) ==  pQuery->substr(0, secondread.length()-pos))
		{
			//return overlapped portion
			mergedseq=secondread.substr(pos);
			return true;
		}
	}

	return false;

}
*/
/*
size_t LongReadExtend::calculateKmerCoverage (const std::string & seq , size_t kmerLength , const BWT* pBWT)
{
	if (seq.length() < kmerLength) return 0;

	size_t cov = 0 ;
	for (size_t i=0; i<=seq.length()-kmerLength;i+=kmerLength/2)
		cov += BWTAlgorithms::countSequenceOccurrences(seq.substr(i,kmerLength) , pBWT );
	
	return cov;
}

// replace each kmer with highest one at each locus
bool LongReadExtend::replaceLowFreqKmer (std::string & seq , size_t kmerLength)
{
	bool changed = false;
	
	for (size_t i=0; i <=seq.length()-kmerLength; i++)
	{
		//Forward kmer should be computed reversely using pRBWT for forward extension
		BWTInterval fwdProbe=BWTAlgorithms::findInterval(indices.pRBWT, reverse(seq.substr(i, kmerLength-1)));
		BWTInterval rvcProbe=BWTAlgorithms::findInterval(indices.pBWT, reverseComplement(seq.substr(i, kmerLength-1)));
		
		size_t maxcov=0;
		for(int j = 1; j < BWT_ALPHABET::size; ++j) //j=A,C,G,T
		{
			char b = BWT_ALPHABET::getChar(j);

			//update forward Interval using extension b
			if(fwdProbe.isValid())
				BWTAlgorithms::updateInterval(fwdProbe, b, indices.pRBWT);

			//update reverse complement Interval using extension rcb
			char rcb=BWT_ALPHABET::getChar(5-i); //T,G,C,A
			if(rvcProbe.isValid())
				BWTAlgorithms::updateInterval(rvcProbe, rcb, indices.pBWT);

			size_t bcount = 0;
			if(fwdProbe.isValid())
				bcount += fwdProbe.size();
			if(rvcProbe.isValid())
				bcount += rvcProbe.size();

			if(bcount > maxcov) {
				maxcov = bcount;
				seq.replace(i+kmerLength-1, 1, 1, b);
				changed = true;
			}
		}
	}
	
	return changed;
}

// Remove leaves with two or more same kmers
void LongReadExtend::removeLeavesByRepeatKmer()
{
	STNodePtrList newLeaves;

	for(STNodePtrList::iterator iter = leaves.begin(); iter != leaves.end(); ++iter)
	{
		std::string STNodeStr = (*iter)->getFullString();
		std::string fwdrepeatunit = STNodeStr.substr(STNodeStr.size()-minOverlap);
		std::string revrepeatunit = reverseComplement(fwdrepeatunit);
		size_t index1=STNodeStr.find(fwdrepeatunit);
		size_t index2=STNodeStr.find(revrepeatunit);

		if(index1 == (STNodeStr.size()- minOverlap) && index2 == std::string::npos)
		{
			newLeaves.push_back(*iter);
		}
	}

	leaves=newLeaves;
}

// Print the string represented by every node
void LongReadExtend::printAll()
{
	std::cout << "Print all: \n";
	pRootNode->printAllStrings("");
}
*/
