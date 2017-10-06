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
//Formula for calculating kmer threshold value ( x,y,z ) = ( isLowCoverage,m_params.PBcoverage,staticKmerSize )
#define FORMULA( x,y,z ) ( (x) ? (0.05776992234 * y - 0.4583043394 * z + 10.19159685) : (0.0710704607 * y - 0.5445663957 * z + 12.26253388) )
//Generic judgment of a kmer ( a,b,c,d ) = ( currentKmerFreqs,dynamicKmerThresholdValue,fwdKmerFreqs,rvcKmerFreqs )
#define OVERTHRESHOLD( a,b,c,d ) ( (a>=b) && (c>=1) && (d>=1) )
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
   
	//Find seeds using fixed or dynamic kmers depending on 1st round or not     
    seedVec = hybridSeedingFromPB(readSeq);
       
	result.Timer_Seed = seedTimer->getElapsedWallTime(); 
    delete seedTimer;
	result.totalSeedNum = seedVec.size();
    

	//Push the first seed into pacbioCorrectedStrs, which will be popped later as source seed
	if(seedVec.size() >= 2 )
	{
		result.correctedLen += seedVec[0].seedStr.length();		
		pacbioCorrectedStrs.push_back(seedVec[0]);		
	}
	else
	{
		//Give up reads with less than 2 seeds
		result.merge = false;
		return result;
	}

    //Rserve sufficient str length for fast append
    pacbioCorrectedStrs.back().seedStr.reserve(readSeq.length());
	
    initCorrect(readSeq, seedVec, pacbioCorrectedStrs, result);

	
	result.merge = true;
	result.totalReadsLen = readSeq.length();
	for(size_t i = 0 ; i < pacbioCorrectedStrs.size() ; i++)
		result.correctedPacbioStrs.push_back(pacbioCorrectedStrs[i].seedStr);
	
	return result;
    
}
/*Legacy code ignored by KuanWeiLee
void PacBioCorrectionProcess::separatebykmer(std::string readid,std::string readSeq,size_t staticKmerSize)
{
    std::string outfilename = "k"+readid + "_" + std::to_string(staticKmerSize) + "mer.sf";
    ofstream outfile (outfilename) ;
    outfile << ">"<<readid << "\n" ;
    for(size_t i = 0 ; i+staticKmerSize <= readSeq.length() ; i++)
	{
		std::string kmer = readSeq.substr(i, staticKmerSize);
    
		size_t fwdKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(kmer, m_params.indices);
		size_t rvcKmerFreqs = BWTAlgorithms::countSequenceOccurrencesSingleStrand(reverseComplement(kmer), m_params.indices);
		size_t currentKmerFreqs = fwdKmerFreqs+rvcKmerFreqs;
        outfile << kmer << "\t" << currentKmerFreqs << "\n" ;
    }
    outfile.close();
    
}
*/
void PacBioCorrectionProcess::initCorrect(std::string& readSeq, std::vector<SeedFeature>& seedVec, std::vector<SeedFeature>& pacbioCorrectedStrs, PacBioCorrectionResult& result)
{


    m_total_FMtime = 0;
    m_total_DPtime = 0;
    FMextendParameters FMextendParameter(m_params.indices,m_params.idmerLength,m_params.maxLeaves,m_params.minKmerLength,m_params.PBcoverage,m_params.ErrorRate,m_params.DebugExtend);
    
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
			target =  seedVec.at(targetSeed + nextTargetSeed);
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
			// if(dis_between_src_target >= (int)m_params.maxSeedInterval&& m_params.PBcoverage >=50) 
				// break;

            
            
			// extension using local kmer hashtable collected from overlapping reads
			std::string mergedseq;
			FMWalkReturnType = extendBetweenSeeds(source, target, readSeq, mergedseq, extendKmerSize, dis_between_src_target,FMextendParameter);
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
                // target =  targetSeed+nextTargetSeed+1<seedVec.size()?seedVec[targetSeed+nextTargetSeed+1]:target;
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


// Dynamic seeding from pacbio reads, v20160517 by Ya.
// identify seeds by dynamic kmers from pacbio reads, 
// which is suitable for the PacBio hybrid error correction,
// where repeat regions require large kmers and error-prone regions require small kmers.

// Seeding by fixed and dynamic kmer size 
std::vector<SeedFeature> PacBioCorrectionProcess::hybridSeedingFromPB(const std::string& readSeq)
{   
    std::vector<SeedFeature> seedVec;
    const size_t staticKmerSize = m_params.kmerLength;
    // Prevention of short reads
	if(readSeq.length() < staticKmerSize) 
		return seedVec;       
    
	//Slide with a fixed kmer. Noted by KuanWeiLee
	std::vector<BWTIntervalPair> FixedMerInterval;
    std::vector<size_t> freqsCount;
    freqsCount.assign(m_params.PBcoverage*2,0);
    for(size_t i = 0 ; i <= readSeq.length() - staticKmerSize ; i++)
	{
		//Collection of fixed kmer interval(s) sliding at every position
		//on the read sequence(index:i [0~(readSeq.length()-staticKmerSize)]). Noted by KuanWeiLee
        std::string kmer = readSeq.substr(i,staticKmerSize);
		BWTInterval fwdInterval = BWTAlgorithms::findInterval(m_params.indices.pRBWT, reverse(kmer));
		BWTInterval rvcInterval = BWTAlgorithms::findInterval(m_params.indices.pBWT, reverseComplement(kmer));
        BWTIntervalPair bip;
		bip.interval[0] = fwdInterval;
        bip.interval[1] = rvcInterval;
        FixedMerInterval.push_back(bip);
		//Statistics of freqsCount(index:currentKmerFreqs)		
		size_t currentKmerFreqs = bip.getFreqs();
		if(currentKmerFreqs < freqsCount.size())
            freqsCount[currentKmerFreqs]++;
		
		//[Debugseed] should be output more dedicately by KuanWeiLee
        if(m_params.DebugSeed)
                    std::cout << i << ": "<< kmer << " total " << currentKmerFreqs <<":" << fwdInterval.size() << ":" << rvcInterval.size() << " <=\n";   
    }       
    //Determine whether the read is of low coverage;
	//Thresholds need further inspectation.Noted & abbreviated by KuanWeiLee
    size_t initKmerThresholdValueWithLowCoverage = FORMULA(true,m_params.PBcoverage,staticKmerSize);
    size_t initKmerThresholdValue = FORMULA(false,m_params.PBcoverage,staticKmerSize);
	bool isLowCoverage = freqsCount[(int)initKmerThresholdValueWithLowCoverage] > freqsCount[(int)initKmerThresholdValue];
	
	//Threshold table (index:j [staticKmerSize~kmerLengthUpperBound]) is calculated by a formula; 
	//there are 2 modes in the formula: on low coverage or not, and the lower bound is 3. Noted by KuanWeiLee
    std::vector<float> kmerThresholdTable;
	const size_t kmerLengthUpperBound = 50;
	kmerThresholdTable.assign(kmerLengthUpperBound+1,0);
	for(size_t j=staticKmerSize ; j <= kmerLengthUpperBound ; j++)
	{        
		float kmerThresholdValue = FORMULA(isLowCoverage,m_params.PBcoverage,j);
		kmerThresholdTable[j] = kmerThresholdValue  < 3 ? 3: kmerThresholdValue;
	}
    
	//[Debugseed] should be output more dedicately.Noted & abbreviated by KuanWeiLee
    if(m_params.DebugSeed)
		std::cout << ( isLowCoverage ? "is low" : "isn't low" ) << std::endl;
    
    size_t lastRepeatSeedPos = 0;
    size_t lastRepeatSeedFreqs = 0;
    //Searching seed; loop through the whole kmers with fixed size on the read sequence. Noted by KuanWeiLee
	/*
	for(size_t i=0; i <= readSeq.length() - staticKmerSize; i++)
	{
		bool isSeed = false;
		char b,rcb;
		std::string kmer = readSeq.substr(i,staticKmerSize);
		size_t dynamicKmerSize = staticKmerSize, dynamicKmerThresholdValue;
		BWTInterval fwdInterval = FixedMerInterval[i].interval[0],
					rvcInterval = FixedMerInterval[i].interval[1];
		size_t fwdKmerFreqs, rvcKmerFreqs, currentKmerFreqs;
		size_t fixedMerFreqs, maxFixedMerFreqs=0;
		size_t seedStartPos = i, seedEndPos;
		float GCratio;
		for(size_t j = i; j <= readSeq.length() - staticKmerSize; j++, isSeed = true)
		{	
			
			if(isSeed)
			{
				b = readSeq[j + staticKmerSize - 1];
				switch(b)
				{
					case 'A' : rcb='T'; break;
					case 'T' : rcb='A'; break;
					case 'C' : rcb='G'; break;
					case 'G' : rcb='C'; break;
				}
				kmer = kmer + b;
				dynamicKmerSize++;
				BWTAlgorithms::updateInterval(fwdInterval,b,m_params.indices.pRBWT);
				BWTAlgorithms::updateInterval(rvcInterval,rcb,m_params.indices.pBWT);
			}
			dynamicKmerThresholdValue = kmerThresholdTable[dynamicKmerSize];
			fwdKmerFreqs = fwdInterval.getFreqs();
			rvcKmerFreqs = rvcInterval.getFreqs();			
			currentKmerFreqs = fwdKmerFreqs + rvcKmerFreqs;
			fixedMerFreqs = FixedMerInterval[j].getFreqs();
			maxFixedMerFreqs = std::max(maxFixedMerFreqs,fixedMerFreqs);
			
			if(!OVERTHRESHOLD(currentKmerFreqs,dynamicKmerThresholdValue,fwdKmerFreqs,rvcKmerFreqs)
				|| isLowComplexity(kmer,GCratio) || dynamicKmerSize > kmerLengthUpperBound)
			{
				if(isSeed)
				{
					dynamicKmerSize--;
					kmer.erase(dynamicKmerSize);
				}
				break;
			}
		}
		if(isSeed)
		{
			bool isRepeat = maxFixedMerFreqs > 5*kmerThresholdTable[staticKmerSize];
			seedEndPos = seedStartPos + dynamicKmerSize - 1;
			SeedFeature newSeed(seedStartPos, kmer, isRepeat, staticKmerSize, m_params.PBcoverage/2);
			newSeed.estimateBestKmerSize(m_params.indices.pBWT);
			newSeed.maxFixedMerFreqs = maxFixedMerFreqs;
			seedVec.push_back(newSeed);
			i = seedEndPos;
		}
	}
	*/
    //Structure needs shrinking. Noted by KuanWeiLee
	for(size_t i = 0 ; i <= readSeq.length() - staticKmerSize ; i++)
    {
        
        std::string kmer = readSeq.substr(i,staticKmerSize);
        size_t dynamicKmerSize = staticKmerSize;
		size_t dynamicKmerThresholdValue = kmerThresholdTable[dynamicKmerSize];
        BWTInterval fwdInterval = FixedMerInterval[i].interval[0];
        BWTInterval rvcInterval = FixedMerInterval[i].interval[1];
		size_t fwdKmerFreqs = fwdInterval.getFreqs();
		size_t rvcKmerFreqs = rvcInterval.getFreqs();
        size_t currentKmerFreqs = fwdKmerFreqs + rvcKmerFreqs;
		size_t fixedMerFreqs = currentKmerFreqs;
                
		//[Debugseed] should be output more dedicately by KuanWeiLee
        if(m_params.DebugSeed)
            std::cout << i << ": " << kmer << "\t" << currentKmerFreqs <<":" << fwdKmerFreqs << ":" << rvcKmerFreqs << "\n"; //debugch
        
		if( OVERTHRESHOLD(currentKmerFreqs,dynamicKmerThresholdValue,fwdKmerFreqs,rvcKmerFreqs) )
		{	
			//Skip low-complexity seeds, e.g., AAAAAAAATAAAA, which are often error seeds
			float GCRatio = 0;
			if(isLowComplexity(kmer, GCRatio))
			{
				bool isPrevSeedCloseToRepeat = !seedVec.empty() && !seedVec.back().isRepeat 
												// remove previous seed if it's too close 
												&& i - seedVec.back().seedEndPos < m_repeat_distance 
												// ensure that the kmer frequency difference is also large
												&& ((float)seedVec.back().endKmerFreq/(float)currentKmerFreqs  < 0.6); 
										
				
				if(isPrevSeedCloseToRepeat)	
					seedVec.pop_back();
				
				continue;
			}
            
            if( isCloseTorepeat(FixedMerInterval,i,m_repeat_distance))
            {
                // std::cout<<  "yyy\n";
                continue;
            }
      
                                   
			size_t seedStartPos = i;
			//The maximun fixed kmer freq would be traced every position. Noted by KuanWeiLee
			size_t maxFixedMerFreqs = fixedMerFreqs;
			
			//Cluster consecutive fixed kmers as one seed as long al possible. Noted by KuanWeiLee
			for(i++ ; i <= readSeq.length() - staticKmerSize ; i++)
			{
                //get next char
                char b = readSeq[i + staticKmerSize -1];
                char rcb;
                switch(b)
                {
                    case 'A': rcb='T'; break;
                    case 'T': rcb='A'; break;
                    case 'C': rcb='G'; break;
                    case 'G': rcb='C'; break;
                }
                
                //Extend kmer length and upadate intervals. Noted by KuanWeiLee
                kmer = kmer + b;
				dynamicKmerSize++;
				dynamicKmerThresholdValue = kmerThresholdTable[dynamicKmerSize];
                BWTAlgorithms::updateInterval(fwdInterval,b,m_params.indices.pRBWT);
                BWTAlgorithms::updateInterval(rvcInterval,rcb,m_params.indices.pBWT);
				fwdKmerFreqs = fwdInterval.getFreqs();
				rvcKmerFreqs = rvcInterval.getFreqs();
                currentKmerFreqs = fwdKmerFreqs + rvcKmerFreqs;				
                fixedMerFreqs = FixedMerInterval[i].getFreqs();				
                
                maxFixedMerFreqs = std::max(maxFixedMerFreqs,fixedMerFreqs);
                
                //[Debugseed] should be output more dedicately by KuanWeiLee
                if(m_params.DebugSeed)
                    std::cout << i << ": "<< kmer << "\t local "<< fixedMerFreqs << " total " << currentKmerFreqs <<":" << fwdKmerFreqs << ":" << rvcKmerFreqs <<" <=\n"; //debugch
                
				//If too long/of low complexity/hitch-hiked/not over seed extension threshold, jump out the [FOOR] loop. Noted by KuanWeiLee				
				if( dynamicKmerSize > kmerLengthUpperBound || isLowComplexity(kmer, GCRatio) || (float)fixedMerFreqs/(float)maxFixedMerFreqs < 0.6 ||
				!(OVERTHRESHOLD(currentKmerFreqs,dynamicKmerThresholdValue,fwdKmerFreqs,rvcKmerFreqs) && (int)fixedMerFreqs >= (int)initKmerThresholdValueWithLowCoverage) )
				{
					kmer.erase(dynamicKmerSize-1);
					dynamicKmerSize--;
					break;
				}				
			}//End clustering consecutive fixed kmers as one seed [FOOR]. Noted by KuanWeiLee	
			size_t seedEndPos = seedStartPos + dynamicKmerSize -1;

			//Upadate seed vector and detect whether the seed is in repeat. Noted by KuanWeiLee
			if(maxFixedMerFreqs > kmerThresholdTable[staticKmerSize]*5)
			{
                
                

                // PB135123_7431.fa TGTAATCAGGCTGAAAA
                bool isPrevSeedBetweenRepeat = seedVec.size()>=2 
										// Previous first seed is not repeat
										&& !seedVec.back().isRepeat 
                                        // but previous second seed and current seed are both repeats
                                        && seedVec[seedVec.size()-2].isRepeat
                                        // and the distance between two repeat seeds are not large
                                        && seedStartPos - seedVec[seedVec.size()-2].seedEndPos < m_repeat_distance*2; 
                                        
                //PB135123_7431.fa: CGGCTTTCTTTAATGAT
                bool isPrevSeedWithinLargeRepeat = seedVec.size()>=3
										// Previous seed is not repeat
										&& !seedVec.back().isRepeat 
                                        // but both previous two seeds are repeats
                                        && seedVec[seedVec.size()-2].isRepeat && seedVec[seedVec.size()-3].isRepeat
                                        // and the distance between two repeat seeds are not large
                                        && seedStartPos - seedVec[seedVec.size()-2].seedEndPos < m_repeat_distance*4; 


                if(isPrevSeedBetweenRepeat || isPrevSeedWithinLargeRepeat)	
                    seedVec.pop_back();

                // PB135123_7431.fa: AAATTTGAAGAGACTCA, CCAGGGTATCTAAATCCTGTTT
                bool isPrevTwoSeedWithinLargeRepeat = seedVec.size()>=4 
										// Previous two seed are not repeat
										&& !seedVec.back().isRepeat && !seedVec[seedVec.size()-2].isRepeat 
                                        // but previous seed is repeats
                                        && seedVec[seedVec.size()-3].isRepeat 
                                        // and the distance between two repeat seeds are not large
                                        && seedStartPos - seedVec[seedVec.size()-3].seedEndPos < m_repeat_distance*4 
                                        && (seedVec.back().seedLength-staticKmerSize<=3 || seedVec[seedVec.size()-2].seedLength-staticKmerSize<=3);

                if(isPrevTwoSeedWithinLargeRepeat)
				{
                    // std::cout << "hahaha" << seedVec.back().isSmall() << "\t" << seedVec.back().seedStr << "\n";                
                    seedVec.pop_back();
                    seedVec.pop_back();
                }
										
				// if(isPrevSeedCloseToRepeat || isPrevSeedBetweenRepeat || isPrevSeedWithinLargeRepeat)	
					// seedVec.pop_back();
				         
                if(maxFixedMerFreqs > 1024)                
					continue;
                
                
                // std::cout<<     "repeats\n";
				i=seedEndPos;
				SeedFeature newSeed(seedStartPos, kmer, true, staticKmerSize, m_params.PBcoverage/2);
				newSeed.estimateBestKmerSize(m_params.indices.pBWT);
                newSeed.maxFixedMerFreqs = maxFixedMerFreqs;
				seedVec.push_back(newSeed);	
				continue;
				//Restart from end of previous repeat seed because multiple repeat seeds may be within the same region (?) Noted by KuanWeiLee	
			}
			else 
			{				
				SeedFeature newSeed(seedStartPos, kmer, false, staticKmerSize, m_params.PBcoverage/2);
				newSeed.estimateBestKmerSize(m_params.indices.pBWT);
                newSeed.maxFixedMerFreqs = maxFixedMerFreqs;
				seedVec.push_back(newSeed);			
			}			
			i = seedEndPos;
		}
	}//End searching seed [FOOR]. Noted by KuanWeiLee
	
	//[Debugseed] should be output more dedicately.Noted & abbreviated by KuanWeiLee
    if(m_params.DebugSeed)
    {
         ofstream outfile ( m_readid+"_seedfile2.out") ;
         // outfile<<">"+m_readid<< "\n";
         for(size_t j=0;j<seedVec.size();j++)
            outfile<< seedVec[j].seedStr << "\t" << seedVec[j].maxFixedMerFreqs<< "\t" << seedVec[j].seedStartPos<< "\n";
         outfile.close();
    }
	return seedVec;
}


//check if in repeat region and there is high variation 
bool PacBioCorrectionProcess::isCloseTorepeat(std::vector<BWTIntervalPair> FixedMerInterval,size_t &currpos,size_t m_repeat_distance)
{
    
    size_t NearbyMaxFreqs = 0;
    //Curr Seed Freqs
    size_t CurrSeedFreqs  = FixedMerInterval[currpos].interval[0].size() + FixedMerInterval[currpos].interval[1].size();
    size_t lowerbound = currpos - m_repeat_distance > 0 ? currpos - m_repeat_distance : 0;
    
    size_t uperbound = currpos + m_repeat_distance > FixedMerInterval.size()-1 ?  FixedMerInterval.size()-1 : currpos + m_repeat_distance;
    
    //Search Maximum Freqs in this region
    for(size_t repeatCheckPos = lowerbound; repeatCheckPos <= uperbound ; repeatCheckPos++ )
    {
        
         size_t CurrPosFreqs = FixedMerInterval[repeatCheckPos].interval[0].size() + FixedMerInterval[repeatCheckPos].interval[1].size();
         //std::cout << CurrPosFreqs << "\n";
         if (CurrPosFreqs > NearbyMaxFreqs) NearbyMaxFreqs = CurrPosFreqs ;
    }
   
    bool isLargeVariationOfFreqs = (float)CurrSeedFreqs/(float)NearbyMaxFreqs < 0.6;
    
    return isLargeVariationOfFreqs;


}

// Perform FMindex extension between source and target seeds
// Return FMWalkReturnType

int PacBioCorrectionProcess::extendBetweenSeeds(SeedFeature& source, SeedFeature& target, std::string& rawSeq, std::string& mergedseq, 
												size_t extendKmerSize, size_t dis_between_src_target,FMextendParameters FMextendParameter)
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
   if(source.isRepeat && ! target.isRepeat)
   {    
        // std::cout<<      "0.0\n";
       
        LongReadSelfCorrectByOverlap OverlapTree(reverseComplement(target.seedStr),reverseComplement(strbetweensrctarget),reverseComplement(srcStr),dis_between_src_target,extendKmerSize,extendKmerSize+2,FMextendParameter,min_SA_threshold);
        FMWalkReturnType = 	OverlapTree.extendOverlap(fmwalkresult);
        fmwalkresult.mergedSeq = reverseComplement(fmwalkresult.mergedSeq);
   
   }
    else
    {
        LongReadSelfCorrectByOverlap OverlapTree(srcStr,strbetweensrctarget,target.seedStr,dis_between_src_target,extendKmerSize,extendKmerSize+2,FMextendParameter,min_SA_threshold);
        FMWalkReturnType = 	OverlapTree.extendOverlap(fmwalkresult);  
      
    }
    
       if(FMWalkReturnType > 0)//extend success by fm extend
        {
            
          mergedseq = fmwalkresult.mergedSeq;
          mergedseq = mergedseq.substr(extendKmerSize);
        }

          m_total_FMtime += FMTimer->getElapsedWallTime() ;
    
     
    delete FMTimer;
    
   
    
    
    
    
    
   
    if(FMWalkReturnType <= 0)
    //v2
    {
        
        
        double identity = 0.66;
         size_t totalMaxFreqs = source.maxFixedMerFreqs + target.maxFixedMerFreqs;
		if(m_params.DebugExtend)
			std::cout << source.maxFixedMerFreqs + target.maxFixedMerFreqs << "	  freqs\n";
		size_t min_call_coverage = 15;
		if (totalMaxFreqs > 50)
		{
			identity =0.7;
			min_call_coverage = totalMaxFreqs*0.4;
		   
		}
		if (totalMaxFreqs > 100)
		{
			identity =0.75;
			min_call_coverage = totalMaxFreqs*0.4;
		   
        }
        Timer* DPTimer = new Timer("DP Time",true);
        std::string rawSubseq = source.seedStr.substr(source.seedStr.length()-extendKmerSize) +  strbetweensrctarget + target.seedStr;
        
        MultipleAlignment maquery = LongReadOverlap::buildMultipleAlignment(rawSubseq,
                                                        extendKmerSize, //m_params.PBKmerLength, 
                                                        extendKmerSize, //m_params.PBKmerLength,
                                                        rawSubseq.length()/10, 
                                                        identity,	// alignment identity < 0.7 are often false positive repeats
                                                        m_params.PBcoverage,
                                                        m_params.indices);
        
        std::string consensus = maquery.calculateBaseConsensus(min_call_coverage, -1);
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
	
	size_t staticKmerSize = m_params.kmerLength;
	
	std::string kmer = readSeq.substr(seedStartPos, staticKmerSize);
	int initKmerFreq = BWTAlgorithms::countSequenceOccurrences(kmer, m_params.indices);
	int prevKmerFreq = initKmerFreq;

	// first kmer is a repeat
	if(initKmerFreq > minRepeatFreq)
	{
		newSeedStartPos = seedStartPos;
		startKmerFreq = initKmerFreq;
	}
	
	
	// identify breakpoints of large freq difference between two kmers	
	for(size_t i=seedStartPos+1 ; i+staticKmerSize-1 <= seedEndPos; i++)
	{
		kmer = readSeq.substr(i, staticKmerSize);
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
				newSeedEndPos = i + staticKmerSize -2;
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

bool PacBioCorrectionProcess::isLowComplexity (const std::string& seq, float& GCratio, float threshold)
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
	GCratio = (float)(countG+countC)/seqLen;
	 return (float) countA/seqLen >= threshold || (float) countT/seqLen >= threshold
			|| (float) countC/seqLen >= threshold || (float) countG/seqLen >= threshold;
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
SeedFeature::SeedFeature(size_t startPos, std::string str, bool repeat, size_t staticKmerSize, size_t repeatCutoff)
	:seedStartPos(startPos), seedStr(str), isRepeat(repeat), minKmerSize(staticKmerSize), freqUpperBound(repeatCutoff),
	freqLowerBound(repeatCutoff/2), stepSize(1)
{
	seedEndPos = seedStartPos + seedStr.length() -1;
	seedLength = seedStr.length();
	startBestKmerSize = endBestKmerSize = staticKmerSize<=seedLength?staticKmerSize:seedLength;
}

// append current seed string with extendedStr
void SeedFeature::append(std::string extendedStr)
{ 
	seedStr += extendedStr;
	seedLength += extendedStr.length();
	seedStartPos += extendedStr.length();
	seedEndPos += extendedStr.length();
}

void SeedFeature::setBestKmerSize(size_t staticKmerSize)
{
	startBestKmerSize = endBestKmerSize = staticKmerSize;
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
