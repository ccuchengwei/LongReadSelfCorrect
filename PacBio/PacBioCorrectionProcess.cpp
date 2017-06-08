///-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioCorrectionProcess.cpp - Self-correction using FM-index walk for PacBio reads
//
#include "PacBioCorrectionProcess.h"
#include "SAIPBSelfCTree.h"
#include "CorrectionThresholds.h"
#include "HashMap.h"
#include <iomanip>
#include <time.h>
#include "SAIPBHybridCTree.h"


#include "LongReadOverlap.h"
#include "Timer.h"


using namespace std;


PacBioCorrectionProcess::PacBioCorrectionProcess(const PacBioCorrectionParameters params) : m_params(params)
{
}

PacBioCorrectionProcess::~PacBioCorrectionProcess()
{

}


// PacBio Self Correction by Ya and YTH, v20151202.
// 1. Identify highly-accurate seeds within PacBio reads
// 2. For each pair of seeds, perform kmer extension using local kmer frequency collected by FM-index extension
PacBioCorrectionResult PacBioCorrectionProcess::PBSelfCorrection(const SequenceWorkItem& workItem)
{	
	// std::cout << workItem.read.id << "\n";
    m_readid =  workItem.read.id;
	PacBioCorrectionResult result;
	
    Timer* seedTimer = new Timer("Seed Time",true);
  

	std::vector<SeedFeature> seedVec, pacbioCorrectedStrs;
	std::string readSeq = workItem.read.seq.toString();

   // separatebykmer(workItem.read.id,readSeq,m_params.kmerLength);
   
	// find seeds using fixed or dynamic kmers depending on 1st round or not
     
    seedVec = hybridSeedingFromPB(readSeq);
    
    
   
     
    
	result.Timer_Seed = seedTimer->getElapsedWallTime(); 
    delete seedTimer;
	result.totalSeedNum = seedVec.size();
    

	// push the first seed into pacbioCorrectedStrs, which will be popped later as source seed
	if(seedVec.size() >= 2)
	{
		result.correctedLen += seedVec.at(0).seedStr.length();
		// if(m_params.isSplit)
			pacbioCorrectedStrs.push_back(seedVec.at(0));
		// else
			// pacbioCorrectedStrs.push_back(readSeq.substr(0, seedVec.at(0).seedEndPos+1));
		
		// if(!m_params.isSplit)
			// pacbioCorrectedStrs.back().seedStr.reserve(readSeq.length()*1.5);
	}
	else
	{
		// give up reads with less than 2 seeds
		result.merge = false;
		return result;
	}
	

    // reserve sufficient str length for fast append
    pacbioCorrectedStrs.back().seedStr.reserve(readSeq.length());
    initCorrect(readSeq, seedVec, pacbioCorrectedStrs, result);

	
	result.merge = true;
	result.totalReadsLen = readSeq.length();
	for(size_t result_count = 0 ; result_count < pacbioCorrectedStrs.size() ; result_count++)
		result.correctedPacbioStrs.push_back(pacbioCorrectedStrs[result_count].seedStr);
	
	return result;
    
}

void PacBioCorrectionProcess::separatebykmer(std::string readid,std::string readSeq,size_t kmerSize)
{
    std::string outfilename = "k"+readid + "_" + std::to_string(kmerSize) + "mer.sf";
    ofstream outfile (outfilename) ;
    outfile << ">"<<readid << "\n" ;
    for(size_t i = 0 ; i+kmerSize <= readSeq.length() ; i++)
	{
		std::string kmer = readSeq.substr(i, kmerSize);
    
		size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
		size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
		size_t kmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
        outfile << kmer << "\t" << kmerFreqs << "\n" ;
    }
    outfile.close();
    
}

void PacBioCorrectionProcess::initCorrect(std::string& readSeq, std::vector<SeedFeature>& seedVec, std::vector<SeedFeature>& pacbioCorrectedStrs, PacBioCorrectionResult& result)
{


    m_total_FMtime = 0;
    m_total_DPtime = 0;
    // FMextendParameters setFMextendParameter();
    
	for(size_t targetSeed = 1 ; targetSeed < seedVec.size() ; targetSeed++)
	{				
		// number of trials of extension to the same target seed
		// size_t numOfTrials = 0;
		
		int FMWalkReturnType = 0;
        // int prevFMWalkReturnType = 0;
		
		// source is increasing because no split in the 1st round, this is slow, better replace with pointer
		SeedFeature source = pacbioCorrectedStrs.back();
		SeedFeature target =  seedVec.at(targetSeed);

		// extension kmer is used for extension using local kmer hashtable collected from overlapping reads
		// default: smaller beset kmer size from both seeds -2
		size_t extendKmerSize = std::min(source.endBestKmerSize, seedVec.at(targetSeed).startBestKmerSize) - 2;
			
            
		// Multiple targets will be tested for FM-index walk from source to target, if target is error seed, until m_params.numOfNextTarget times.
		for(int nextTargetSeed = 0 ; nextTargetSeed < (int)m_params.numOfNextTarget && targetSeed + nextTargetSeed < seedVec.size() ; nextTargetSeed++)
		{
			// std::cout << "======= " << result.totalWalkNum << " =======\n";
			
			// extendKmerSize should be large for repeat seeds
			if((source.isRepeat || target.isRepeat) )
			{
				extendKmerSize = std::min((int)source.seedLength, (int)target.seedLength);
				if(int(extendKmerSize) > m_params.kmerLength+2) 
						extendKmerSize = m_params.kmerLength+2;
			}

			// Estimate distance between source and target, but this may over-estimate due to more insertion errors
			// Note that source seed has been updated and no long stands for the original seed, which is seedVec[targetSeed-1]
			int dis_between_src_target = target.seedStartPos - seedVec.at(targetSeed-1).seedStartPos - seedVec.at(targetSeed-1).seedStr.length();
		
			// skip seeds with large distance in between for speedup
			if(dis_between_src_target >= (int)m_params.maxSeedInterval&& m_params.PBcoverage >=50) 
				break;



			// extension using local kmer hashtable collected from overlapping reads
			std::string mergedseq;
			FMWalkReturnType = extendBetweenSeeds(source, target, readSeq, mergedseq, extendKmerSize, dis_between_src_target);
            if(m_params.DebugExtend)
            std::cout<< targetSeed << " \t" <<FMWalkReturnType<< " result of extension\n"<< mergedseq << "\n"; //debugch
            // outfile << targetSeed << " \t" <<FMWalkReturnType << "\n";
			if(FMWalkReturnType > 0)
			{
				// FMWalk success
				// size_t extendStartPos = source.seedLength;
				// std::string extendedStr = mergedseq.substr(extendStartPos);
				
				// append extended string into last corrected seed string and update related seed attributes
				// pacbioCorrectedStrs.back().append(extendedStr);
				pacbioCorrectedStrs.back().append(mergedseq);
				// the last seed will become new source and should be updated
				pacbioCorrectedStrs.back().endBestKmerSize = target.endBestKmerSize;
				pacbioCorrectedStrs.back().isRepeat = target.isRepeat;

				// result statistics
				// result.correctedLen += extendedStr.length();
				result.correctedLen += mergedseq.length();
                if(FMWalkReturnType ==1)
                    result.correctedNum++;
                else
                    result.DPNum++;
				result.seedDis += dis_between_src_target;
				
				// jump to nextTargetSeed+1 if more than one target was tried and succeeded
				targetSeed = targetSeed + nextTargetSeed;
                result.totalWalkNum++;
				break;
			}
            /*
			else
			{
				// return <0: give up this source seed
				int ActionFlag = FMWalkFailedActions(extendKmerSize, numOfTrials, source, target, FMWalkReturnType, prevFMWalkReturnType);
				if(ActionFlag <0)
					break;
				// return 0: retry the same target
				else if(ActionFlag == 0)
					nextTargetSeed--;
				// return >0: move on to next target
				else
				{
					// target =  targetSeed+nextTargetSeed+1<seedVec.size()?seedVec[targetSeed+nextTargetSeed+1]:target;
					// extendKmerSize = std::min(source.endBestKmerSize, target.startBestKmerSize) - 2;
				}
			}
           */

			// prevFMWalkReturnType = FMWalkReturnType;
		}// end of next target seed
		
		// All targets failure: 
		// 0: seed inter-distance too large
		// -1: kmer extension failed at later stage close to target
		// -4: kmer extension failed at early stage close to source
		// -2: exceed depth
		// -3: exceed leaves
		if(FMWalkReturnType <= 0)
		{
			// push seedVec[targetSeed] into results, which will become new source in the next iteration
			result.seedDis += seedVec[targetSeed].seedStartPos - seedVec[targetSeed-1].seedStartPos - seedVec[targetSeed-1].seedStr.length();
			result.correctedLen += seedVec[targetSeed].seedStr.length();

			// retain uncorrected parts of reads
			if(!m_params.isSplit)
			{
				size_t startPos = seedVec[targetSeed-1].seedStartPos + seedVec[targetSeed-1].seedStr.length();
				size_t extendedLen = seedVec[targetSeed].seedStartPos + seedVec[targetSeed].seedStr.length() - startPos;

				pacbioCorrectedStrs.back().append(readSeq.substr(startPos,extendedLen));
				pacbioCorrectedStrs.back().endBestKmerSize = seedVec.at(targetSeed).endBestKmerSize;
				pacbioCorrectedStrs.back().isRepeat = seedVec.at(targetSeed).isRepeat;
			}
			else
			{
				// split original read into seeds and discard uncorrected parts of reads
				pacbioCorrectedStrs.push_back(seedVec[targetSeed]);
			}
			
			// statistics of FM extension
			if(FMWalkReturnType == -1 || FMWalkReturnType == -4)
				result.highErrorNum++;
			else if(FMWalkReturnType == -2)
				result.exceedDepthNum++;
			else if(FMWalkReturnType == -3)
				result.exceedLeaveNum++;
            if(FMWalkReturnType != 0)
            result.totalWalkNum++;
		}
		
		
	}// end of each target seed
    
    result.Timer_FM = m_total_FMtime;
    result.Timer_DP = m_total_DPtime;
    // cout<< m_total_FMtime<<endl;
    // outfile.close();
            
            // std::cout<< alnscore.first << " | " << alnscore.second << "kk\n";
    
    
}
/*
FMextendParameters PacBioCorrectionProcess:setFMextendParameter()
{
    
    FMextendParameters FMextendParametersFile;
    if(m_params.PBcoverage < 60)
        FMextendParametersFile.
    
    
    
    
}*/
void PacBioCorrectionProcess::realCorrect(std::string& readSeq, std::vector<SeedFeature>& seedVec, std::vector<SeedFeature>& pacbioCorrectedStrs, PacBioCorrectionResult& result)
{
	// for each pair of seeds, perform kmer extension using local kmer frequency
	for(size_t targetSeed = 1 ; targetSeed < seedVec.size() ; targetSeed++)
	{				
		// number of trials of extension to the same target seed
		size_t numOfTrials = 0;
		
		int FMWalkReturnType = 0, prevFMWalkReturnType = 0;
		SeedFeature source = pacbioCorrectedStrs.back();
		
		//PB36993_4517, GTAATAAGGAATATCTCAATTTT is a repeat seed with large frequency leading to -2
		// re-estimate the best kmer size if source becomes larger
		if(source.endKmerFreq>90)
			source.estimateBestKmerSize(m_params.indices.pBWT);
		
		SeedFeature target =  seedVec.at(targetSeed);

		// small kmer is used for extension using kmer hashtable from overlapping reads
		size_t extendKmerSize = std::min(source.endBestKmerSize, seedVec.at(targetSeed).startBestKmerSize) - 2;

		// Multiple targets will be tested for FM-index walk from source to target, if target is error seed, until m_params.numOfNextTarget times.
		for(int nextTargetSeed = 0 ; nextTargetSeed < (int)m_params.numOfNextTarget && targetSeed + nextTargetSeed < seedVec.size() ; nextTargetSeed++)
		{
			// std::cout << "======= " << result.totalWalkNum << " =======\n";
			
			// Estimate distance between source and target, but this may over-estimate due to more insertion errors
			// Note that source seed has been updated and no long stands for the original seed, which is seedVec[targetSeed-1]
			int dis_between_src_target = target.seedStartPos - seedVec.at(targetSeed-1).seedStartPos - seedVec.at(targetSeed-1).seedStr.length();
					
			// skip seeds with large distance in between for speedup
			// if(dis_between_src_target >= (int)m_params.maxSeedInterval) 
				// break;

			// PB159615_16774.fa require smaller kmer
			// else if(m_params.isFirst && dis_between_src_target >= 50 && !source.isRepeat && !target.isRepeat && extendKmerSize>11)
				// extendKmerSize -= 2;
		

			std::string mergedseq;
			FMWalkReturnType = extendBetweenSeeds(source, target, readSeq, mergedseq, extendKmerSize, dis_between_src_target);

			if(FMWalkReturnType > 0)
			{
				// FMWalk success
				// size_t extendStartPos = source.seedLength;
				// std::string extendedStr = mergedseq.substr(extendStartPos);
				
				// append extended string into last corrected seed string and update related seed attributes
				// pacbioCorrectedStrs.back().append(extendedStr);
				pacbioCorrectedStrs.back().append(mergedseq);
				
				// the last seed will become new source and should be updated
				pacbioCorrectedStrs.back().endBestKmerSize = target.endBestKmerSize;
				pacbioCorrectedStrs.back().isRepeat = target.isRepeat;

				// result statistics
				// result.correctedLen += extendedStr.length();
				result.correctedLen += mergedseq.length();
				result.correctedNum++;
				result.seedDis += dis_between_src_target;
				
				// jump to nextTargetSeed+1 if more than one target was tried and succeeded
				targetSeed = targetSeed + nextTargetSeed;
				break;
			}
			else
			{ 
				// return <0: give up this source seed
				int ActionFlag = FMWalkFailedActions(extendKmerSize, numOfTrials, source, target, FMWalkReturnType, prevFMWalkReturnType);
				if(ActionFlag <0)
					break;
				// return 0: retry the same target
				else if(ActionFlag == 0)
					nextTargetSeed--;
				// return >0: move on to next target
				else
				{
					// target =  targetSeed+nextTargetSeed+1<seedVec.size()?seedVec[targetSeed+nextTargetSeed+1]:target;
					// extendKmerSize = std::min(source.endBestKmerSize, target.startBestKmerSize) - 2;
				}
			}
			
			prevFMWalkReturnType = FMWalkReturnType;
		}// end of next target seed
		
		// All targets failure: 
		// 0: seed inter-distance too large
		// -1: kmer extension failed at later stage close to target
		// -4: kmer extension failed at early stage close to source
		// -2: exceed depth
		// -3: exceed leaves
		if(FMWalkReturnType <= 0)
		{
			// push seedVec[targetSeed] into results, which will become new source in the next iteration
			result.seedDis += seedVec[targetSeed].seedStartPos - seedVec[targetSeed-1].seedStartPos - seedVec[targetSeed-1].seedStr.length();
			result.correctedLen += seedVec[targetSeed].seedStr.length();

			// retain uncorrected parts of reads
			if(!m_params.isSplit)
			{
				size_t startPos = seedVec[targetSeed-1].seedStartPos + seedVec[targetSeed-1].seedStr.length();
				size_t extendedLen = seedVec[targetSeed].seedStartPos + seedVec[targetSeed].seedStr.length() - startPos;

				pacbioCorrectedStrs.back().append(readSeq.substr(startPos,extendedLen));
				pacbioCorrectedStrs.back().endBestKmerSize = seedVec.at(targetSeed).endBestKmerSize;
				pacbioCorrectedStrs.back().isRepeat = seedVec.at(targetSeed).isRepeat;
			}
			else
			{
				// split original read into seeds and discard uncorrected parts of reads
				pacbioCorrectedStrs.push_back(seedVec[targetSeed]);
			}
			
			// statistics of FM extension
			if(FMWalkReturnType == -1 || FMWalkReturnType == -4)
				result.highErrorNum++;
			else if(FMWalkReturnType == -2)
				result.exceedDepthNum++;
			else if(FMWalkReturnType == -3)
				result.exceedLeaveNum++;
			
		}
		
		//result.totalWalkNum++;
	}// end of each target seed
}


int PacBioCorrectionProcess::checkseedcorrect(std::vector<SeedFeature> seeds,std::string currseedStr,size_t currseedStartpos)
{
	int *correctnum = new int[currseedStr.length()];
	for(size_t i=0;i<currseedStr.length();i++)correctnum[i]=0;
	size_t runtime=1;
	size_t correctSeedlength=0;
   
	if(seeds.size()>=2&&!seeds.at(seeds.size()-1).isRepeat) runtime=2;
	const size_t kmerSize = m_params.kmerLength;
	for(;runtime>0;runtime--)
	{
	
	std::string preseedStr = seeds.at(seeds.size()-runtime).seedStr;
	std::string initKmer = preseedStr.substr(preseedStr.length() - kmerSize);
	BWTInterval fwdInterval=BWTAlgorithms::findInterval(m_params.indices.pRBWT, reverse(initKmer));
	BWTInterval rvcInterval=BWTAlgorithms::findInterval(m_params.indices.pBWT, reverseComplement(initKmer));
	size_t extendlength = currseedStartpos-seeds.at(seeds.size()-runtime).seedStartPos+1+seeds.at(seeds.size()-runtime).seedStr.length();
	            // cout<<"preseed:"<<preseedStr<<"\ncurrseedStr:"<<currseedStr<<"\n======================\n";
				for(int64_t fwdRootIndex = fwdInterval.lower;fwdInterval.isValid()&& fwdRootIndex <= fwdInterval.upper ;fwdRootIndex++)
				  {
					 std::string currentFwdKmer = initKmer;
					 std::string extendkmer = currentFwdKmer ;
					 int64_t fwdIndex = fwdRootIndex;
					 correctSeedlength=0;
					 for(int64_t currentLength = (int64_t)kmerSize; currentLength <= 1.5*extendlength; currentLength++)
					 {
						 
						 char b = m_params.indices.pRBWT->getChar(fwdIndex);
			             if(b == '$'||kmerSize+correctSeedlength>currseedStr.length()) break;
						 currentFwdKmer = currentFwdKmer.substr(1) + b ;
						 extendkmer = extendkmer+b;
						 if(currentFwdKmer == currseedStr.substr(correctSeedlength,kmerSize))
						 {
							 correctnum[correctSeedlength]++;
							 correctSeedlength++;
							 // cout << extendkmer <<"\n"; 
						 }
						 fwdIndex = m_params.indices.pRBWT->getPC(b) + m_params.indices.pRBWT->getOcc(b, fwdIndex - 1);
					 }
					 // cout << extendkmer <<"\n";
				  }
				  // cout<<"===============\nRvccurrseed:"<<reverseComplement(currseedStr.substr(0,kmerSize))<<"\n";
				 
				  for(int64_t rvcRootIndex=rvcInterval.lower;rvcRootIndex <= rvcInterval.upper && rvcInterval.isValid();rvcRootIndex++)
				{
					std::string currentRvcKmer = reverseComplement(initKmer);
					std::string extendkmer = currentRvcKmer ;
					int64_t rvcIndex = rvcRootIndex;
					 correctSeedlength=0;
					for(int64_t currentLength = (int64_t)kmerSize; currentLength <= 1.5*extendlength;currentLength++)
					{
					   char b = m_params.indices.pBWT->getChar(rvcIndex);
			           if(b == '$'||kmerSize+correctSeedlength>currseedStr.length()) break;
					   currentRvcKmer = b+currentRvcKmer.substr(0,kmerSize-1) ;
					   
					   if( currentRvcKmer == reverseComplement(currseedStr.substr(correctSeedlength,kmerSize)))
					   {
						   correctnum[correctSeedlength]++;
							correctSeedlength++;
					   }
						rvcIndex = m_params.indices.pBWT->getPC(b) + m_params.indices.pBWT->getOcc(b, rvcIndex - 1);
					}
					// cout << extendkmer <<"<=\n";
				}
	}	
	 int correctThreshold = 3;

				for(correctSeedlength=0;kmerSize+correctSeedlength<=currseedStr.length();correctSeedlength++)
				{
					if(correctnum[correctSeedlength]<correctThreshold&&correctSeedlength>0)
					{
						
						delete [] correctnum;
						return kmerSize+correctSeedlength-1 ;
					}
					else if (correctnum[correctSeedlength]<correctThreshold&&correctSeedlength==0)
					{
						delete [] correctnum;
						return -1;
					}
				}
				
				delete [] correctnum;
				return kmerSize+correctSeedlength-1;
}


// dynamic seeding from pacbio reads, v20160517 by Ya.
// identify seeds by dynamic kmers from pacbio reads, 
// which is suitable for the PacBio hybrid error correction,
// where repeat regions require large kmers and error-prone regions require small kmers.

// seeding by fixed and dynamic kmer size 
std::vector<SeedFeature> PacBioCorrectionProcess::hybridSeedingFromPB(const std::string& readSeq)
{   
    std::vector<SeedFeature> seedVec;
    const size_t kmerSize = m_params.kmerLength;
    // prevention of short reads
	if(readSeq.length() < kmerSize) return seedVec;

    
    std::vector<BWTIntervalPair> FixedMerInterval; 
    
    
    size_t kmerThresholdValueWithLowCoverage = 0.05776992234 * m_params.PBcoverage - 0.4583043394 * kmerSize + 10.19159685;
    size_t kmerThresholdValue = 0.0710704607 * m_params.PBcoverage - 0.5445663957 * kmerSize + 12.26253388;
     
    std::vector<size_t> freqscount;
    
    freqscount.resize(m_params.PBcoverage*2,0);
    for(size_t i = 0 ; i+kmerSize <= readSeq.length() ; i++)
	{
        std::string kmer = readSeq.substr(i,kmerSize);
		BWTInterval fwdInterval = BWTAlgorithms::findInterval(m_params.indices.pRBWT, reverse(kmer));
		BWTInterval rvcInterval = BWTAlgorithms::findInterval(m_params.indices.pBWT, reverseComplement(kmer));
        BWTIntervalPair bip;
        bip.interval[0] = fwdInterval;
        bip.interval[1] = rvcInterval;
        
		size_t kmerFreqs = (fwdInterval.isValid()?fwdInterval.size():0) + (rvcInterval.isValid()?rvcInterval.size():0);
        
        FixedMerInterval.push_back(bip);
        if(kmerFreqs < freqscount.size())
            freqscount.at(kmerFreqs) += 1;
    
    }
    
    
    bool isLowCoverage = false;
    //define low coverage read
    if(freqscount.at(kmerThresholdValueWithLowCoverage) > freqscount.at(kmerThresholdValue))
        isLowCoverage= true;
 
    
    std::vector<float> kmerThreshold;
	kmerThreshold.resize(51,0);
	for(size_t fixedmerSize=kmerSize ; fixedmerSize<51 ; fixedmerSize++)
	{
        float kmerThresholdValue;
        if(!isLowCoverage)
            // kmerThresholdValue=(-0.5*fixedmerSize+16.17)*((float)m_params.PBcoverage/60);//threshold of normal coverage
            kmerThresholdValue = 0.0710704607*m_params.PBcoverage-0.5445663957*kmerSize+ 12.26253388;
        else    
            // kmerThresholdValue = (-0.43*fixedmerSize+14.1)*((float)m_params.PBcoverage/60);//threshold of low coverage
            kmerThresholdValue =  0.05776992234*m_params.PBcoverage-0.4583043394*kmerSize+10.19159685;
		kmerThreshold.at(fixedmerSize) = kmerThresholdValue  < 3 ? 3: kmerThresholdValue;
	}
    
    //debug mode
    if(m_params.DebugSeed)
    {
        if(isLowCoverage)
            std::cout<< "is low\n";
        else
            std::cout<< "isn't low\n";
    }
    
    
    //start seqrching seed  
    for(size_t i = 0 ; i + kmerSize <= readSeq.length() ; i++)
    {
        //find initial kmer frequency
        std::string kmer = readSeq.substr(i,kmerSize);
        size_t dynamicKmerSize = kmerSize;
        BWTInterval fwdInterval = FixedMerInterval.at(i).interval[0];
        BWTInterval rvcInterval = FixedMerInterval.at(i).interval[1];
        size_t fwdKmerFreqs = fwdInterval.isValid()?fwdInterval.size():0;
        size_t rvcKmerFreqs = rvcInterval.isValid()?rvcInterval.size():0;
        size_t kmerFreqs = fwdKmerFreqs + rvcKmerFreqs;
        size_t dynamicKmerThreshold = kmerThreshold.at(dynamicKmerSize);
        
        //debug info
        if(m_params.DebugSeed)
            std::cout << i << ": " << kmer << "\t" << kmerFreqs <<":" << fwdKmerFreqs << ":" << rvcKmerFreqs << "\n"; //debugch
        
        
		if(kmerFreqs >= dynamicKmerThreshold && fwdKmerFreqs>=1 && rvcKmerFreqs>=1)
		{	
			// skip low-complexity seeds, e.g., AAAAAAAATAAAA, which are often error seeds
			float GCRatio = 0;
			if(isLowComplexity(kmer, GCRatio))
			{
				bool isPrevSeedCloseToRepeat = !seedVec.empty() && !seedVec.back().isRepeat 
										// remove previous seed if it's too close to the current repeat within kmerSize
										&& i - seedVec.back().seedEndPos < m_repeat_distance 
										// ensure that the kmer frequency difference is also large
										&& seedVec.back().seedLength-kmerSize<=3; 
				
				if(isPrevSeedCloseToRepeat)	
					seedVec.pop_back();
				
				continue;
			}
			size_t seedStartPos = i, 
			seedLen = 0;
			
			// Group consecutive solid kmers into one seed if possible
			size_t maxKmerFreq = kmerFreqs;
			for(i++ ; i + kmerSize <= readSeq.length(); i++)
			{
                //get next char
                char b = readSeq.at( i + kmerSize -1);
                char rcb;
                switch(b)
                {
                    case 'A': rcb='T'; break;
                    case 'T': rcb='A'; break;
                    case 'C': rcb='G'; break;
                    case 'G': rcb='C'; break;
                }
                
                //update interval
                kmer = kmer + b;
                BWTAlgorithms::updateInterval(fwdInterval,b,m_params.indices.pRBWT);
                BWTAlgorithms::updateInterval(rvcInterval,rcb,m_params.indices.pBWT);
                kmerFreqs = (fwdInterval.isValid()?fwdInterval.size():0) + (rvcInterval.isValid()?rvcInterval.size():0);
                fwdKmerFreqs = fwdInterval.isValid()?fwdInterval.size():0;
                rvcKmerFreqs = rvcInterval.isValid()?rvcInterval.size():0;
                
                
                size_t FixedMerFreqs = FixedMerInterval.at(i).interval[0].size() + FixedMerInterval.at(i).interval[1].size();
                
                
                dynamicKmerSize++;
                
                
                // assert(dynamicKmerSize<=kmerThreshold.size());
                //if seed is too long ,cut off
                if(dynamicKmerSize >= kmerThreshold.size()) break;
                //skip low-complexity
                if( isLowComplexity(kmer, GCRatio) ) break;
                
                
                //set threshold with corresponding length
                dynamicKmerThreshold = kmerThreshold.at(dynamicKmerSize);
                maxKmerFreq = std::max(maxKmerFreq,FixedMerFreqs);
                
                
                if(m_params.DebugSeed)
                    std::cout << i << ": "<< kmer << "\t local "<< FixedMerFreqs << " total " << kmerFreqs <<":" << fwdKmerFreqs << ":" << rvcKmerFreqs << " || " << dynamicKmerThreshold <<" <=\n"; //debugch
                
                 
                if(((int)kmerFreqs >= (int)dynamicKmerThreshold) && (int)fwdKmerFreqs>=1 && (int)rvcKmerFreqs>=1 && (int)FixedMerFreqs >= (int)kmerThresholdValueWithLowCoverage)
                    seedLen++;
                //seed error
                else
                {

                    dynamicKmerSize--;
                    dynamicKmerThreshold = kmerThreshold.at(dynamicKmerSize);
      
                    break;
                }
                
                /*
                if(kmerFreqs >= kmerThreshold.at(kmerSize) && FixedMerInterval.at(i).interval[0].size()>=1 && FixedMerInterval.at(i).interval[1].size()>=1)
                    seedLen++;

                else
                
                    break;
                
                */
                 
                
			} //second for end
			
			size_t seedEndPos = seedStartPos + seedLen + kmerSize-1;
			// skip contaminated seeds
			


			

			// this is a repeat seed
			if(maxKmerFreq > kmerThreshold.at(kmerSize)*5)	//repeat seed
			{
                
                // For seeds within repeats, error seeds will still exceed the cutoff, e.g., 12 11 15 60 65 70 20 19..
                // refine the exact boundary by finding the segments with highest kmer frequency 60 65 70
                // std::pair<size_t, size_t> kmerFreqPair = refineRepeatSeed(readSeq, seedStartPos, seedEndPos,kmerThreshold.at(kmerSize));

                // PB135123_7431.fa, error seeds close to repeats are also with high frequency
                bool isPrevSeedCloseToRepeat = !seedVec.empty() && !seedVec.back().isRepeat 
                                        // remove previous seed if it's too close 
                                        && seedStartPos - seedVec.back().seedEndPos < m_repeat_distance 
                                        // ensure that the kmer frequency difference is large
                                        && (std::abs((int)seedVec.back().endKmerFreq-(int)maxKmerFreq)>kmerThreshold.at(kmerSize)*3); 
                
                if(isPrevSeedCloseToRepeat)	
                    seedVec.pop_back();

                // PB135123_7431.fa TGTAATCAGGCTGAAAA
                bool isPrevSeedBetweenRepeat = seedVec.size()>=2 && !seedVec.back().isRepeat // previous seed is not repeat
                                        // but previous previous seed and current seed are both repeats
                                        && seedVec.at(seedVec.size()-2).isRepeat
                                        // and the distance between two repeat seeds are not large
                                        && seedStartPos - seedVec.at(seedVec.size()-2).seedEndPos < m_repeat_distance*2; 
                                        
                //PB135123_7431.fa: CGGCTTTCTTTAATGAT
                bool isPrevSeedWithinLargeRepeat = seedVec.size()>=3 && !seedVec.back().isRepeat // previous seed is not repeat
                                        // but both previous two seeds are repeats
                                        && seedVec.at(seedVec.size()-2).isRepeat && seedVec.at(seedVec.size()-3).isRepeat
                                        // and the distance between two repeat seeds are not large
                                        && seedStartPos - seedVec.at(seedVec.size()-2).seedEndPos < m_repeat_distance*4; 


                if(isPrevSeedBetweenRepeat || isPrevSeedWithinLargeRepeat)	
                    seedVec.pop_back();

                // PB135123_7431.fa: AAATTTGAAGAGACTCA, CCAGGGTATCTAAATCCTGTTT
                bool isPrevTwoSeedWithinLargeRepeat = seedVec.size()>=4 && !seedVec.back().isRepeat && !seedVec.at(seedVec.size()-2).isRepeat // previous two seed are not repeat
                                        // but previous seed is repeats
                                        && seedVec.at(seedVec.size()-3).isRepeat 
                                        // and the distance between two repeat seeds are not large
                                        && seedStartPos - seedVec.at(seedVec.size()-3).seedEndPos < m_repeat_distance*4 
                                        && (seedVec.back().seedLength-kmerSize<=3 || seedVec.at(seedVec.size()-2).seedLength-kmerSize<=3);

                if(isPrevTwoSeedWithinLargeRepeat){
                    // std::cout << "hahaha" << seedVec.back().isSmall() << "\t" << seedVec.back().seedStr << "\n";
                
                    seedVec.pop_back();
                    seedVec.pop_back();
                }
										
				// if(isPrevSeedCloseToRepeat || isPrevSeedBetweenRepeat || isPrevSeedWithinLargeRepeat)	
					// seedVec.pop_back();
				
				// Unsolved cases
				// PB38933_7552.fa
				// 1414: GTTCAGCGGAAATTTTC 1:1:0
				// 1415: TTCAGCGGAAATTTTCC 20:11:9
				// 1416: TCAGCGGAAATTTTCCA 1:1:0
				// Seed:   17      TTCAGCGGAAATTTTCC
				
				/*
				if(!seedVec.empty())
				{
					int checkresult=checkseedcorrect(seedVec,readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1),i);
					// cout<<seedEndPos-seedStartPos+1<<"\t"<<checkresult<<"\t"<<i<<endl;
					
					 if(checkresult>=17)
					 {   
						 seedEndPos=seedStartPos + checkresult -1;
						 
						}
					 if(checkresult==-1){ 
							
							continue;
						}
				}*/
                
                if(maxKmerFreq > kmerThreshold.at(kmerSize) * 11)
                
				continue;
                
                
                
				i=seedEndPos;
				SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), true, kmerSize, m_params.PBcoverage/2);
				newSeed.estimateBestKmerSize(m_params.indices.pBWT);
                newSeed.maxFixedMerFreqs = maxKmerFreq;
				seedVec.push_back(newSeed);	
				continue;
				
							
			}
			else 
			{
				bool isCloseToPrevRepeatSeed = !seedVec.empty() && seedVec.back().isRepeat 
												&& (seedStartPos - seedVec.back().seedEndPos <= kmerSize);
									
									
				// most error seeds are not too large or with lower frequency
				// PB135123_7431.fa: AAAACTTCGCAGTGAAC is not error but discarded
				// && seedEndPos+1-seedStartPos-kmerSize < 7;
												
				// if(seedVec.empty() || !seedVec.back().isRepeat || (seedStartPos - seedVec.back().seedEndPos > kmerSize) )
				if(seedVec.empty() || !isCloseToPrevRepeatSeed)
                   
				{
					// push concatenated seeds into seed vector
					SeedFeature newSeed(seedStartPos, readSeq.substr(seedStartPos, seedEndPos-seedStartPos+1), false, kmerSize, m_params.PBcoverage/2);
					// newSeed.setBestKmerSize(kmerSize);
					newSeed.estimateBestKmerSize(m_params.indices.pBWT);
                    newSeed.maxFixedMerFreqs = maxKmerFreq;
					seedVec.push_back(newSeed);
				}
			

			}
			
			
			// restart from end of previous repeat seed because multiple repeat seeds may be within the same region
			i=seedEndPos;

		}// end of sufficient kmerThreshold

	}// end of for
    // if(m_params.DebugSeed)
    // {
         // ofstream outfile ( m_readid+"_seedfile2.out") ;
         // outfile<<">"+m_readid<< "\n";
         // for(size_t j=0;j<seedVec.size();j++)

            // outfile<< seedVec.at(j).seedStr << "\t" << seedVec.at(j).maxFixedMerFreqs<< "\t" << seedVec.at(j).seedStartPos<< "\n";
            


         // outfile.close();
    // }
	return seedVec;
}



// Perform FMindex extension between source and target seeds
// Return FMWalkReturnType

int PacBioCorrectionProcess::extendBetweenSeeds(SeedFeature& source, SeedFeature& target, std::string& rawSeq, std::string& mergedseq, 
												size_t extendKmerSize, size_t dis_between_src_target)
{
   // size_t srcKmerSize = std::max(source.endBestKmerSize, extendKmerSize);
   std::string srcStr = source.seedStr.substr(source.seedStr.length()-extendKmerSize);
   std::string strbetweensrctarget = rawSeq.substr(target.seedStartPos-dis_between_src_target,dis_between_src_target);

   int min_SA_threshold =3;
    //v1
    if (m_params.PBcoverage > 60) min_SA_threshold = int(m_params.PBcoverage/60)*3;

    
    
    FMWalkResult2 fmwalkresult;
    Timer* FMTimer = new Timer("FM Time",true);
    int FMWalkReturnType =0;

        LongReadSelfCorrectByOverlap OverlapTree(source.seedStr,strbetweensrctarget,target.seedStr,dis_between_src_target,extendKmerSize,extendKmerSize-2,extendKmerSize+2,m_params.indices,m_params.PBcoverage,m_params.maxLeaves,min_SA_threshold,m_params.DebugExtend);
        FMWalkReturnType = 	OverlapTree.extendOverlap(fmwalkresult);
        if(FMWalkReturnType > 0)//extend success by fm extend
        {
            
            mergedseq = fmwalkresult.mergedSeq;

            if(!mergedseq.empty())
                mergedseq = mergedseq.substr(extendKmerSize);
        }
          m_total_FMtime += FMTimer->getElapsedWallTime() ;
    
     
    delete FMTimer;
    
   
    if(FMWalkReturnType <= 0)
    //v2
    {
        Timer* DPTimer = new Timer("DP Time",true);
        std::string rawSubseq = source.seedStr.substr(source.seedStr.length()-extendKmerSize) +  strbetweensrctarget + target.seedStr;
        
        MultipleAlignment maquery = LongReadOverlap::buildMultipleAlignment(rawSubseq,
                                                        extendKmerSize, //m_params.PBKmerLength, 
                                                        extendKmerSize, //m_params.PBKmerLength,
                                                        rawSubseq.length()/10, 
                                                        0.7,	// alignment identity < 0.7 are often false positive repeats
                                                        m_params.PBcoverage,
                                                        m_params.indices);
        
        std::string consensus = maquery.calculateBaseConsensus(200, -1);
        // std::cout << rawSubseq << "   raw\n";
        // std::cout << ">" << consensus.length() <<"\n" << consensus << "   <-- consensus"<<endl;
        if(maquery.getNumRows() <= 3)
                return FMWalkReturnType;
        mergedseq = consensus;
        if(!mergedseq.empty())
                mergedseq = mergedseq.substr(extendKmerSize);
        
        m_total_DPtime += DPTimer->getElapsedWallTime() ;
        delete DPTimer;
        
        FMWalkReturnType = 2;
    
    }
    
    
	return FMWalkReturnType;
}
// refine seed interval using larger kmer
std::pair<size_t, size_t> PacBioCorrectionProcess::refineRepeatSeed(const std::string readSeq, size_t& seedStartPos, size_t& seedEndPos,size_t normal_freqs)
{
	// initially set to max unisnged int value
	size_t newSeedStartPos = (size_t)-1;
	size_t newSeedEndPos = (size_t)-1;
	size_t startKmerFreq=0, endKmerFreq=0;
	
	const int minRepeatFreq = normal_freqs*3, minFreqDiff = 30;
	
	size_t kmerSize = m_params.kmerLength;
	
	std::string kmer = readSeq.substr(seedStartPos, kmerSize);
	int initKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
	int prevKmerFreq = initKmerFreq;

	// first kmer is a repeat
	if(initKmerFreq > minRepeatFreq)
	{
		newSeedStartPos = seedStartPos;
		startKmerFreq = initKmerFreq;
	}
	
	
	// identify breakpoints of large freq difference between two kmers	
	for(size_t i=seedStartPos+1 ; i+kmerSize-1 <= seedEndPos; i++)
	{
		kmer = readSeq.substr(i, kmerSize);
		int currKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);

		// std::cout << i << ": " << kmer << "\t" << currKmerFreq << "\n";

		// error kmers within repeats often lead to large freq diff
		bool isLargeFreqDiff = currKmerFreq - prevKmerFreq > minFreqDiff;
		
		// PB36993_4517.fa, TTATGTAAGGAGTATTTGAT
		// some error kmers are with moderate frequency and no freq diff can be observed
		// pick up the first repeat kmer as starting point
		bool isRepeatKmer = (newSeedStartPos == (size_t)-1) && (currKmerFreq >= (int)minRepeatFreq);
		if(isLargeFreqDiff || isRepeatKmer)
		{
			// capture the highest repeat start pos
			bool isBetterRepeatKmer = (startKmerFreq!=0 && currKmerFreq > (int)startKmerFreq);
			if(newSeedStartPos == (size_t)-1 || isBetterRepeatKmer)
			{
				newSeedStartPos = i;
				startKmerFreq = currKmerFreq;
			}
		}
			
		// repeat end is reached
		// PB36993_4517.fa, AGGCTTGTCTGTAATCGGG
		if(prevKmerFreq - currKmerFreq > minFreqDiff /*|| currKmerFreq < minFreqDiff*/)
		{
			// do not enter unless start pos was found
			// if(newSeedStartPos != (size_t)-1)
			// {
				newSeedEndPos = i + kmerSize -2;
				endKmerFreq = prevKmerFreq;
				break;
			// }
		}
			
		prevKmerFreq = currKmerFreq;
	}
	
	if(newSeedStartPos == (size_t)-1)
	{
		newSeedStartPos = seedStartPos;
		startKmerFreq = initKmerFreq;
	}
		
	if(newSeedEndPos == (size_t)-1)
	{
		newSeedEndPos = seedEndPos;
		endKmerFreq = prevKmerFreq;
	}
	
	// std::cout << newSeedStartPos << "\t" << newSeedEndPos << "\n";

	seedStartPos = newSeedStartPos;
	seedEndPos = newSeedEndPos;
	return std::make_pair(startKmerFreq, endKmerFreq);
}

// return <0: give up and break
// return 0: retry the same target
// return >0: continue to next target
int PacBioCorrectionProcess::FMWalkFailedActions(size_t& extendKmerSize, size_t& numOfTrials, 
								SeedFeature& source, SeedFeature& target, int FMWalkReturnType, int prevFMWalkReturnType)
{
	numOfTrials++;
	// extension failed due to insufficient kmers, reduce large and small kmer sizes
	if(FMWalkReturnType==-1 || FMWalkReturnType==-4)
	{
		// kmers have been enlarged due to repeats, shrink will lead to infinite loop
		if(prevFMWalkReturnType==-3 )
			return -1;
		
		// PB36993_4517.fa, AGGCTTGTCTGTAATCGGG
		if(m_params.isFirst /*&& (source.isRepeat || target.isRepeat)*/)
			return -1;		
		
		source.endBestKmerSize -=2;
		
		target.startBestKmerSize -=2;
		
		extendKmerSize -= 2;

		// std::cout << source.endBestKmerSize << "\t" << target.startBestKmerSize << "\n";

		
		// don't aggressively reduce kmer in the 1st found where most kmers are errors
		if(m_params.isFirst && (source.endBestKmerSize < 15 || target.startBestKmerSize < 15))
			return -1;
			
		if(source.endBestKmerSize < 11 || target.startBestKmerSize < 11 || extendKmerSize < 9)
			return -1;
			
		return 0;
	}
	
	// increase extendKmerSize for reducing repeats
	else if(FMWalkReturnType==-3)
	{
		if(prevFMWalkReturnType==-4 || prevFMWalkReturnType==-1)
			return -1;

		// exponential growth is required in super large repeats. Otherwise speed is too slow
		source.endBestKmerSize += pow(2, numOfTrials+1);
		target.startBestKmerSize += pow(2, numOfTrials+1);
		extendKmerSize += pow(2, numOfTrials+1);
				
		// bug: PB7017_14581_0_14647.fa
		// extendKmerSize is less than seedLength , dunno why
		if(source.seedLength < source.endBestKmerSize || target.seedLength < target.startBestKmerSize ||
			source.seedLength < extendKmerSize || target.seedLength < extendKmerSize )
			return -1;
		
		return 0;
	}
	else if(FMWalkReturnType==-2)
	{
		// probably chimera, need more observations
		// largeKmerSize = m_params.kmerLength;
		// extendKmerSize = largeKmerSize - 2;
		return 1;
	}
	
	return 1;
}

bool PacBioCorrectionProcess::isLowComplexity (std::string seq, float & GCratio, float threshold)
{
	size_t seqLen = seq.length();
	size_t countG =0 ;
	size_t countC =0 ;
	size_t countT =0 ;
	size_t countA =0 ;

	for (size_t i=0; i<seqLen; i++)
	{
		switch(seq[i]){
			case 'A': countA ++ ;break;
			case 'T': countT ++ ;break;
			case 'C': countC ++ ;break;
			case 'G': countG ++ ;break;
			default:  assert(false);
		}
	}

	GCratio = (float)(countG+countC)/seqLen ;

	if (  ((float) countA/seqLen >= threshold ) || ((float) countT/seqLen >= threshold)
			|| ((float) countC/seqLen >= threshold ) || ((float) countG/seqLen >= threshold) )
		return true;

	return false;

}





//
//
//
PacBioCorrectionPostProcess::PacBioCorrectionPostProcess(std::ostream* pCorrectedWriter,
std::ostream* pDiscardWriter,
const PacBioCorrectionParameters params) :
m_pCorrectedWriter(pCorrectedWriter),
m_pDiscardWriter(pDiscardWriter),
m_params(params),
m_totalReadsLen(0),
m_correctedLen(0),
m_totalSeedNum(0),
m_totalWalkNum(0),
m_correctedNum(0),
m_highErrorNum(0),
m_exceedDepthNum(0),
m_exceedLeaveNum(0),
m_DPNum(0),
m_seedDis(0),
m_Timer_Seed(0),
m_Timer_FM(0),
m_Timer_DP(0)
{
}

//
PacBioCorrectionPostProcess::~PacBioCorrectionPostProcess()
{
	if(m_totalWalkNum>0 && m_totalReadsLen>0)
	{
		std::cout << std::endl;
		std::cout << "totalReadsLen: " << m_totalReadsLen << ", ";
		std::cout << "correctedLen: " << m_correctedLen << ", ratio: " << (float)(m_correctedLen)/m_totalReadsLen << "." << std::endl;
		std::cout << "totalSeedNum: " << m_totalSeedNum << "." << std::endl;
		std::cout << "totalWalkNum: " << m_totalWalkNum << ", ";
		std::cout << "FMNum: " << m_correctedNum << ", ratio: " << (float)(m_correctedNum*100)/m_totalWalkNum << "%." << std::endl;
		std::cout << "DPNum: " << m_DPNum << ", ratio: " << (float)(m_DPNum*100)/m_totalWalkNum << "%." << std::endl;
        std::cout << "highErrorNum: " << m_highErrorNum << ", ratio: " << (float)(m_highErrorNum*100)/m_totalWalkNum << "%." << std::endl;
		std::cout << "exceedDepthNum: " << m_exceedDepthNum << ", ratio: " << (float)(m_exceedDepthNum*100)/m_totalWalkNum << "%." << std::endl;
		std::cout << "exceedLeaveNum: " << m_exceedLeaveNum << ", ratio: " << (float)(m_exceedLeaveNum*100)/m_totalWalkNum << "%." << std::endl;
		std::cout << "disBetweenSeeds: " << m_seedDis/m_totalWalkNum << std::endl << std::endl;
        std::cout << "Time of searching Seeds: " << m_Timer_Seed  << std::endl;
        std::cout << "Time of searching FM: " << m_Timer_FM  << std::endl;
        std::cout << "Time of searching DP: " << m_Timer_DP  << std::endl;
	}
}


// Writting results for kmerize and validate
void PacBioCorrectionPostProcess::process(const SequenceWorkItem& item, const PacBioCorrectionResult& result)
{	
	if (result.merge)
	{
		m_totalReadsLen += result.totalReadsLen;
		m_correctedLen += result.correctedLen;
		m_totalSeedNum += result.totalSeedNum;
		m_totalWalkNum += result.totalWalkNum;
		m_correctedNum += result.correctedNum;
		m_highErrorNum += result.highErrorNum;
		m_exceedDepthNum += result.exceedDepthNum;
		m_exceedLeaveNum += result.exceedLeaveNum;
        m_DPNum += result.DPNum;
		m_seedDis += result.seedDis;
        m_Timer_Seed += result.Timer_Seed;
        m_Timer_FM += result.Timer_FM;
        m_Timer_DP += result.Timer_DP;
		//cout << result.correctSequence.toString();
		/*SeqItem mergeRecord;
		stringstream ss;
		ss << item.read.id << "_before_len:" << result.correctSequence.toString().length();
		mergeRecord.id = ss.str();
		mergeRecord.seq = result.correctSequence;
		mergeRecord.write(*m_pCorrectedWriter);*/
		
		for(size_t i = 0 ; i < result.correctedPacbioStrs.size() ; i++)
		{
			SeqItem mergeRecord2;
			std::stringstream ss2;
			ss2 << item.read.id << "_" << i << "_" << result.correctedPacbioStrs[i].toString().length();
			mergeRecord2.id = ss2.str();
			mergeRecord2.seq = result.correctedPacbioStrs[i];
			mergeRecord2.write(*m_pCorrectedWriter);
		}
	}
	else
	{
		// write into discard.fa
		SeqItem mergeRecord2;
		mergeRecord2.id = item.read.id;
		mergeRecord2.seq = item.read.seq;
		mergeRecord2.write(*m_pDiscardWriter);
	}
}

/***************************/
/*** Seed Feature Body *****/
/***************************/
SeedFeature::SeedFeature(size_t startPos, std::string str, bool repeat, size_t kmerSize, size_t repeatCutoff)
	:seedStartPos(startPos), seedStr(str), isRepeat(repeat), freqUpperBound(repeatCutoff), freqLowerBound(repeatCutoff/2), 
	minKmerSize(kmerSize), stepSize(1)
{
	seedEndPos = seedStartPos + seedStr.length() -1;
	seedLength = seedStr.length();
	startBestKmerSize = endBestKmerSize = kmerSize<=seedLength?kmerSize:seedLength;
}

// append current seed string with extendedStr
void SeedFeature::append(std::string extendedStr)
{ 
	seedStr += extendedStr;
	seedLength += extendedStr.length();
	seedStartPos += extendedStr.length();
	seedEndPos += extendedStr.length();
}

void SeedFeature::setBestKmerSize(size_t kmerSize)
{
	startBestKmerSize = endBestKmerSize = kmerSize;
}

void SeedFeature::estimateBestKmerSize(const BWT* pBWT)
{			
	std::string kmerStr = seedStr.substr(0, startBestKmerSize);
	startKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);

	if(startKmerFreq > freqUpperBound)
		increaseStartKmerSize(pBWT);
	else if(startKmerFreq < freqLowerBound)
		decreaseStartKmerSize(pBWT);
		
	kmerStr = seedStr.substr(seedLength-endBestKmerSize);
	endKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);

	if(endKmerFreq > freqUpperBound)
		increaseEndKmerSize(pBWT);
	else if(endKmerFreq < freqLowerBound)
		decreaseEndKmerSize(pBWT);
	
}
	
//estimate kmer size
void SeedFeature::increaseStartKmerSize(const BWT* pBWT)
{
	while(startKmerFreq > freqUpperBound && startBestKmerSize <= seedLength - stepSize)
	{
		startBestKmerSize+=stepSize;
		std::string kmerStr = seedStr.substr(0, startBestKmerSize);
		startKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}
	
	// over increase kmer size
	if(startKmerFreq < freqLowerBound)
	{
		startBestKmerSize-=stepSize;
		std::string kmerStr = seedStr.substr(0, startBestKmerSize);
		startKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}
}

void SeedFeature::decreaseStartKmerSize(const BWT* pBWT)
{
	while(startKmerFreq < freqLowerBound && startBestKmerSize > minKmerSize)
	{
		startBestKmerSize-=stepSize;
		std::string kmerStr = seedStr.substr(0, startBestKmerSize);
		startKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}

	// over reduce kmer size
	if(startKmerFreq>freqUpperBound)
	{
		startBestKmerSize+=stepSize;
		std::string kmerStr = seedStr.substr(0, startBestKmerSize);
		startKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}
}

//estimate kmer size
void SeedFeature::increaseEndKmerSize(const BWT* pBWT)
{
	while(endKmerFreq > freqUpperBound && endBestKmerSize <= seedLength - stepSize)
	{
		endBestKmerSize+=stepSize;
		assert(seedLength >= endBestKmerSize);
		std::string kmerStr = seedStr.substr(seedLength - endBestKmerSize);
		endKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}
	
	if(endKmerFreq < freqLowerBound)
	{
		endBestKmerSize-=stepSize;
		std::string kmerStr = seedStr.substr(seedLength - endBestKmerSize);
		endKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}
}

void SeedFeature::decreaseEndKmerSize(const BWT* pBWT)
{
	while(endKmerFreq < freqLowerBound && endBestKmerSize > minKmerSize)
	{
		endBestKmerSize -= stepSize;
		std::string kmerStr = seedStr.substr(seedLength - endBestKmerSize);
		endKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}
	
	if(endKmerFreq > freqUpperBound)
	{
		endBestKmerSize += stepSize;
		std::string kmerStr = seedStr.substr(seedLength - endBestKmerSize);
		endKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmerStr, pBWT);
	}
}
