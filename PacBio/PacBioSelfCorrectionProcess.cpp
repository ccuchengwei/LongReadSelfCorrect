///-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioSelfCorrectionProcess.cpp - Self-correction using FM-index walk for PacBio reads
//
#include "PacBioSelfCorrectionProcess.h"
#include "SAIPBSelfCTree.h"
#include "CorrectionThresholds.h"
#include "HashMap.h"
#include <iomanip>
#include <time.h>
#include "SAIPBHybridCTree.h"
#include "LongReadOverlap.h"
#include "Timer.h"

//#include "KmerDistribution.h"
//Formula for calculating kmer threshold value ( x,y,z ) = ( isLowCoverage,m_params.PBcoverage,staticKmerSize )
#define FORMULA( x,y,z ) ( (x) ? (0.05776992234f * y - 0.4583043394f * z + 10.19159685f) : (0.0710704607f * y - 0.5445663957f * z + 12.26253388f) )
#define FORMULA2( x,y ) ( 8.641517857 * x - 54.54571078 * y + 1197.001707 )
//Generic judgment of a kmer ( a,b,c,d ) = ( currentKmerFreqs,dynamicKmerThresholdValue,fwdKmerFreqs,rvcKmerFreqs )
#define OVERTHRESHOLD( a,b,c,d ) ( (a>=b) && (c>=1) && (d>=1) )

// PacBio Self Correction by Ya and YTH, v20151202.
// 1. Identify highly-accurate seeds within PacBio reads
// 2. For each pair of seeds, perform kmer extension using local kmer frequency collected by FM-index extension
PacBioSelfCorrectionResult PacBioSelfCorrectionProcess::process(const SequenceWorkItem& workItem)
{	
	PacBioSelfCorrectionResult result;
    result.readid = workItem.read.id;
	result.kd.setID(result.readid);
	
	std::string readSeq = workItem.read.seq.toString();	
	std::vector<SeedFeature> seedVec = hybridSeedingFromPB(readSeq, result), pacbioCorrectedStrs;
	//Push the first seed into pacbioCorrectedStrs, which will be popped later as source seed
	if	(
		seedVec.size() >= 2
		&& !m_params.OnlySeed
		)
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
	for(std::vector<SeedFeature>::iterator iter = pacbioCorrectedStrs.begin(); iter != pacbioCorrectedStrs.end(); iter++)
		result.correctedPacbioStrs.push_back((*iter).seedStr);
	
	return result;
    
}
void PacBioSelfCorrectionProcess::initCorrect(std::string& readSeq, std::vector<SeedFeature>& seedVec, std::vector<SeedFeature>& pacbioCorrectedStrs, PacBioSelfCorrectionResult& result)
{
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
			FMWalkReturnType = extendBetweenSeeds(source, target, readSeq, mergedseq, extendKmerSize, dis_between_src_target, FMextendParameter, result);
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
                pacbioCorrectedStrs.back().maxFixedMerFreqs = target.maxFixedMerFreqs;
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

    // cout<< m_total_FMtime<<endl;
    // outfile.close();
	// std::cout<< alnscore.first << " | " << alnscore.second << "kk\n";
    
    
}


// Seeding by fixed and dynamic kmer size. Noted by KuanWeiLee 20171027
std::vector<SeedFeature> PacBioSelfCorrectionProcess::hybridSeedingFromPB(const std::string& readSeq, PacBioSelfCorrectionResult& result)
{
    Timer* seedTimer = new Timer("Seed Time",true);
    std::vector<SeedFeature> initSeedVec,finalSeedVec;
	std::set<size_t> hitchhikingSeedSet;				//finalSeedVec == initSeedVec - hitchhikingSeedSet
    const size_t staticKmerSize = m_params.kmerLength;
	if(readSeq.length() < staticKmerSize) 
		return initSeedVec;       
    
	//Part 1 : Slide with a fixed kmer. Noted by KuanWeiLee
	std::vector<BistrandBWTInterval> FixedMerInterval;
    std::vector<size_t> freqsCount;
    freqsCount.assign(m_params.PBcoverage*2,0);
	
    for(size_t i = 0 ; i <= readSeq.length() - staticKmerSize ; i++)
	{
		//(1) : Collection of fixed kmer interval(s) at every position on the read sequence
		//(index:i [0~(readSeq.length()-staticKmerSize)]). Noted by KuanWeiLee
        std::string kmer = readSeq.substr(i,staticKmerSize);
		BistrandBWTInterval object;
		object.fwdInterval = BWTAlgorithms::findInterval(m_params.indices.pRBWT, reverse(kmer));
		object.rvcInterval = BWTAlgorithms::findInterval(m_params.indices.pBWT, reverseComplement(kmer));
		FixedMerInterval.push_back(object);
		
		//(2) : Statistics of freqsCount(index:currentKmerFreqs)		
		size_t currentKmerFreqs = object.getFreqs();
		//size_t currentKmerFreqs = bip.getFreqs();
		if(currentKmerFreqs < freqsCount.size())
            freqsCount[currentKmerFreqs]++;
		result.kd.add(currentKmerFreqs);
		//[Debugseed] should be output more dedicately by KuanWeiLee
		/*
		if(m_params.DebugSeed)
		{
			std::cout << i << ": "<< kmer << " total " << currentKmerFreqs << ":";
			std::cout << fwdInterval.getFreqs() << ":" << rvcInterval.getFreqs() << " <=\n";
		}
		*/
    }       
    //Determine whether the read is of low coverage;
	//Thresholds need further inspectation.Noted by KuanWeiLee
	//Threshold table mode of low coverage is canceled. Noted by KuanWeiLee 20171028
    bool isLowCoverage = false;
	{
		//float initKmerThresholdValueWithLowCoverage = FORMULA(true,m_params.PBcoverage,staticKmerSize);
		//float initKmerThresholdValue = FORMULA(false,m_params.PBcoverage,staticKmerSize);
		float initKmerThresholdValue = FORMULA2(m_params.PBcoverage, staticKmerSize);
		//result.kd.computeKDAttributes(1.0f);
		result.kd.computeKDAttributes(initKmerThresholdValue);
		//isLowCoverage = freqsCount[(int)initKmerThresholdValueWithLowCoverage] > freqsCount[(int)initKmerThresholdValue];
	}
	//Part 2 : Get threshold table (index:k [staticKmerSize~kmerLengthUpperBound]); 
	//there are 2 modes in the formula: of low coverage or not, and the lower bound is 5. Noted by KuanWeiLee
    std::vector<float> kmerThresholdTable;
	const int kmerLengthUpperBound = 50;
	kmerThresholdTable.assign(kmerLengthUpperBound+1,0);
	//float weight = 1.f;
	float weight = m_params.PBcoverage / 7.5f * result.kd.getQuartile(1);
	weight = (weight == 0.f ? 1.f : sqrt(weight));
	for(int k = staticKmerSize; k <= kmerLengthUpperBound; k++)
	{        
		//float kmerThresholdValue = FORMULA(isLowCoverage,m_params.PBcoverage,k) * weight;
		float kmerThresholdValue = FORMULA2(m_params.PBcoverage ,k);
		kmerThresholdTable[k] = kmerThresholdValue  < 5 ? 5 : kmerThresholdValue;
	}
	//Part 3 : Search seeds; slide through the read sequence with dynamic kmers. Noted by KuanWeiLee
	for(size_t startPos = 0; startPos <= readSeq.length() - staticKmerSize; startPos++)
	{
		bool isSeed = false;
		std::string kmer = readSeq.substr(startPos,staticKmerSize);
		size_t dynamicKmerSize = staticKmerSize;
		BWTInterval fwdInterval = FixedMerInterval[startPos].fwdInterval,
					rvcInterval = FixedMerInterval[startPos].rvcInterval;
		size_t fwdKmerFreqs, rvcKmerFreqs, currentKmerFreqs;
		size_t fixedMerFreqs, maxFixedMerFreqs = FixedMerInterval[startPos].getFreqs(); 
		float initKmerThresholdValue = kmerThresholdTable[staticKmerSize], dynamicKmerThresholdValue;
		size_t seedStartPos = startPos, seedEndPos;
		float GCratio;
		for(size_t movePos = startPos; movePos <= readSeq.length() - staticKmerSize; movePos++, isSeed = true)
		{	
			
			if(isSeed)
			{
				char b = readSeq[movePos + staticKmerSize - 1];
				char rcb = complement(b);
				kmer += b;
				dynamicKmerSize++;
				BWTAlgorithms::updateInterval(fwdInterval,b,m_params.indices.pRBWT);
				BWTAlgorithms::updateInterval(rvcInterval,rcb,m_params.indices.pBWT);
			}
			dynamicKmerThresholdValue = kmerThresholdTable[dynamicKmerSize];
			fwdKmerFreqs = fwdInterval.getFreqs();
			rvcKmerFreqs = rvcInterval.getFreqs();			
			currentKmerFreqs = fwdKmerFreqs + rvcKmerFreqs;
			fixedMerFreqs = FixedMerInterval[movePos].getFreqs();
			//[Debugseed] should be output more dedicately by KuanWeiLee
			/*
			if(m_params.DebugSeed)
			{
				if (isSeed)
				{
					std::cout << movePos << ": "<< kmer << "\t local "<< fixedMerFreqs << " total " << currentKmerFreqs << ":";
					std::cout << fwdKmerFreqs << ":" << rvcKmerFreqs << " <=" << std::endl;
				}
				else
				{
					std::cout << startPos << ": " << kmer << "\t" << currentKmerFreqs << ":";
					std::cout << fwdKmerFreqs << ":" << rvcKmerFreqs << std::endl;
				}
			}
			*/
			//Gerneral seed extension strategy. Noted by KuanWeiLee
			if  (
				   dynamicKmerSize > kmerLengthUpperBound													//1.over length				
				|| isLowComplexity(kmer,GCratio)															//2.low complexity
				|| fixedMerFreqs < initKmerThresholdValue													//3.fixed kmer frequency
				|| !OVERTHRESHOLD(currentKmerFreqs,dynamicKmerThresholdValue,fwdKmerFreqs,rvcKmerFreqs)		//4.dynamic kmer frequency
				)
			{
				if(isSeed)
				{
					kmer.erase(--dynamicKmerSize);
					seedEndPos = seedStartPos + dynamicKmerSize - 1;
					startPos = seedEndPos;
				}
				break;
			}
			//The order between general seed extentsion and kmer hitchhike strategy may make difference but need advanced observation. Noted by KuanWeiLee 20171027
			//Kmer Hitchhike strategy.
			if	(
				maxFixedMerFreqs > 5*initKmerThresholdValue 
				&& (float)fixedMerFreqs / (float)maxFixedMerFreqs < m_params.hhRatio
				)																							//5.hitchhiking kmer (HIGH-->LOW)
			{
				kmer.erase(--dynamicKmerSize);
				seedEndPos = seedStartPos + dynamicKmerSize - 1;
				startPos = seedEndPos + 1;
				break;
			}
			else if	(
					fixedMerFreqs > 5*initKmerThresholdValue
					&& (float)fixedMerFreqs / (float)maxFixedMerFreqs > m_params.r_hhRatio
					)																						//6.hitchhiking kmer (LOW-->HIGH)
			{
				isSeed = false;
				startPos = movePos - 1;
				break;
			}
			
			maxFixedMerFreqs = std::max(maxFixedMerFreqs,fixedMerFreqs);
		}
		if(isSeed)
		{
			bool isRepeat = maxFixedMerFreqs > 5*initKmerThresholdValue;
			
			//Seed Hitchhike part 1 : mark hitchhiking seeds
			std::vector<SeedFeature>::reverse_iterator iterPrevSeed = initSeedVec.rbegin();
			bool isHitchhiked = false;
			while(iterPrevSeed != initSeedVec.rend() && seedStartPos - (*iterPrevSeed).seedEndPos < m_params.repaetDistance)
			{
				if(!isHitchhiked)
					isHitchhiked = (*iterPrevSeed).isRepeat && (float)(*iterPrevSeed).maxFixedMerFreqs / (float)maxFixedMerFreqs > m_params.r_hhRatio;
				bool prevIsHitchhiked = isRepeat && (float)(*iterPrevSeed).maxFixedMerFreqs / (float)maxFixedMerFreqs < m_params.hhRatio;
				size_t index = initSeedVec.size() -  (iterPrevSeed - initSeedVec.rbegin() + 1);
				if(prevIsHitchhiked)
					hitchhikingSeedSet.insert(index);
				iterPrevSeed++;
			}
			if (isHitchhiked)
				continue;
			
			SeedFeature newSeed(seedStartPos, kmer, isRepeat, staticKmerSize, m_params.PBcoverage/2);
			newSeed.estimateBestKmerSize(m_params.indices.pBWT);
			newSeed.maxFixedMerFreqs = maxFixedMerFreqs;
			initSeedVec.push_back(newSeed);
		}
	}
	//Seed Hitchhike part 2 : delete hitchhiking seeds
	finalSeedVec.reserve(initSeedVec.size());
	{
		size_t index = 0;
		for(std::vector<SeedFeature>::iterator iterSeed = initSeedVec.begin(); iterSeed != initSeedVec.end(); iterSeed++, index++)
			if(hitchhikingSeedSet.count(index) == 0)
				finalSeedVec.push_back((*iterSeed));
	}
	//[Debugseed] Output seeds & kd statistics for each reads.Noted by KuanWeiLee
    if(m_params.DebugSeed)
    {
        std::ostream* pSeedWriter = createWriter(m_params.directory + "seed/" + result.readid + ".seed");
		for(std::vector<SeedFeature>::iterator iterSeed = finalSeedVec.begin(); iterSeed != finalSeedVec.end(); iterSeed++)
			*pSeedWriter << (*iterSeed).seedStr << "\t" << (*iterSeed).maxFixedMerFreqs << "\t" 
			<< (*iterSeed).seedStartPos << "\t" << ((*iterSeed).isRepeat ? "Yes" : "No") << "\n";
		delete pSeedWriter;
		
		std::ostream* pKdWriter = createWriter(m_params.directory + "kmer-stat", std::ios_base::app);
		result.kd.write(*pKdWriter,KmerDistribution::TYPE::ATTRIBUTE);
		delete pKdWriter;
	}
	
	result.totalSeedNum = finalSeedVec.size();
	result.Timer_Seed = seedTimer->getElapsedWallTime(); 
    delete seedTimer;
	return finalSeedVec;
}


//check if in repeat region and there is high variation 
//isCloseTorepeat
/*
bool PacBioSelfCorrectionProcess::isCloseTorepeat(std::vector<BWTIntervalPair> FixedMerInterval,size_t &currpos)
{
    
    size_t NearbyMaxFreqs = 0;
    //Curr Seed Freqs
    size_t CurrSeedFreqs  = FixedMerInterval[currpos].interval[0].size() + FixedMerInterval[currpos].interval[1].size();
    size_t lowerbound = currpos - m_params.repaetDistance > 0 ? currpos - m_params.repaetDistance : 0;
    
    size_t uperbound = currpos + m_params.repaetDistance > FixedMerInterval.size()-1 ?  FixedMerInterval.size()-1 : currpos + m_params.repaetDistance;
    
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
*/

// Perform FMindex extension between source and target seeds
// Return FMWalkReturnType

int PacBioSelfCorrectionProcess::extendBetweenSeeds(SeedFeature& source, SeedFeature& target, std::string& rawSeq, std::string& mergedseq, 
												size_t extendKmerSize, size_t dis_between_src_target,FMextendParameters FMextendParameter, PacBioSelfCorrectionResult& result)
{
   // size_t srcKmerSize = std::max(source.endBestKmerSize, extendKmerSize);
   std::string srcStr = source.seedStr.substr(source.seedStr.length()-extendKmerSize);
   // std::cout<<    source.seedStr << "\n" <<  target.seedStr <<" 0.0\n";
   std::string strbetweensrctarget = rawSeq.substr(target.seedStartPos-dis_between_src_target,dis_between_src_target);

   int min_SA_threshold =3;
    //v1
    if (m_params.PBcoverage > 60) min_SA_threshold = int(m_params.PBcoverage/60)*3;

    
    
    FMWalkResult2 fmwalkresult;
    Timer* FMTimer = new Timer("FM Time",true);
    int FMWalkReturnType =0;
     // std::cout<<  source.maxFixedMerFreqs  <<"||" <<target.maxFixedMerFreqs    <<" test\n";
   if(source.isRepeat && source.maxFixedMerFreqs > target.maxFixedMerFreqs  )
   {    
         // if(extendKmerSize > target.seedStr.length()) extendKmerSize = target.seedStr.length();
        LongReadSelfCorrectByOverlap OverlapTree(reverseComplement(target.seedStr),reverseComplement(strbetweensrctarget),reverseComplement(srcStr),dis_between_src_target,extendKmerSize,extendKmerSize+2,FMextendParameter,min_SA_threshold);
     
       FMWalkReturnType = OverlapTree.extendOverlap(fmwalkresult);
        fmwalkresult.mergedSeq = reverseComplement(fmwalkresult.mergedSeq) + target.seedStr.substr(extendKmerSize);
  
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

    result.Timer_FM += FMTimer->getElapsedWallTime() ;
    
     
    delete FMTimer;
    
   
    
    
    
    
   
   
    if	(
		FMWalkReturnType <= 0
		&& !m_params.NoDp
		)
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
        
        result.Timer_DP += DPTimer->getElapsedWallTime() ;
        delete DPTimer;
        
        FMWalkReturnType = 2;
    
    }
    
     
	return FMWalkReturnType;
}
// refine seed interval using larger kmer
//refineRepeatSeed
/*
std::pair<size_t, size_t> PacBioSelfCorrectionProcess::refineRepeatSeed(const std::string readSeq, size_t& seedStartPos, size_t& seedEndPos,size_t normal_freqs)
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
		//if(prevKmerFreq - currKmerFreq > minFreqDiff || currKmerFreq < minFreqDiff)
		if(prevKmerFreq - currKmerFreq > minFreqDiff)
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
*/
// return <0: give up and break
// return 0: retry the same target
// return >0: continue to next target

//FMWalkFailedActions
/*
int PacBioSelfCorrectionProcess::FMWalkFailedActions(size_t& extendKmerSize, size_t& numOfTrials, 
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
		if(m_params.isFirst)
		//if(m_params.isFirst && (source.isRepeat || target.isRepeat))
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
*/
bool PacBioSelfCorrectionProcess::isLowComplexity (const std::string& seq, float& GCratio, float threshold)
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
PacBioSelfCorrectionPostProcess::PacBioSelfCorrectionPostProcess(std::string correctFile,
std::string discardFile,
const PacBioSelfCorrectionParameters params):
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
	m_pCorrectWriter = createWriter(correctFile);
	m_pDiscardWriter = createWriter(discardFile);
	if(m_params.DebugSeed)
		m_pKdWriter = createWriter(m_params.directory + "kmer-dist");
}

//
PacBioSelfCorrectionPostProcess::~PacBioSelfCorrectionPostProcess()
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
	if(m_params.DebugSeed)
	{
		m_kd.write(*m_pKdWriter);
		delete m_pKdWriter;
	}
	delete m_pCorrectWriter;
	delete m_pDiscardWriter;
}


// Writting results for kmerize and validate
void PacBioSelfCorrectionPostProcess::process(const SequenceWorkItem& item, const PacBioSelfCorrectionResult& result)
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
		mergeRecord.write(*m_pCorrectWriter);*/
		
		for(size_t i = 0 ; i < result.correctedPacbioStrs.size() ; i++)
		{
			SeqItem mergeRecord2;
			std::stringstream ss2;
			ss2 << item.read.id << "_" << i << "_" << result.correctedPacbioStrs[i].toString().length();
			mergeRecord2.id = ss2.str();
			mergeRecord2.seq = result.correctedPacbioStrs[i];
			mergeRecord2.write(*m_pCorrectWriter);
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
	m_kd += result.kd;
}