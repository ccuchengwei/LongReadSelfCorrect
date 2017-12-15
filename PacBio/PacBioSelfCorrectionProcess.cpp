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
#include "KmerThresholdTable.h"

//Generic judgment of a kmer ( a,b,c,d ) = ( currentKmerFreqs,dynamicKmerThresholdValue,fwdKmerFreqs,rvcKmerFreqs )
#define OVERTHRESHOLD( a,b,c,d ) ( (a>=b) && (c>=1) && (d>=1) )
#define MIN(a,b) ( a<b ? a : b)
#define MAX(a,b) ( a>b ? a : b)

// PacBio Self Correction by Ya and YTH, v20151202.
// 1. Identify highly-accurate seeds within PacBio reads
// 2. For each pair of seeds, perform kmer extension using local kmer frequency collected by FM-index extension
PacBioSelfCorrectionResult PacBioSelfCorrectionProcess::process(const SequenceWorkItem& workItem)
{	
	PacBioSelfCorrectionResult result;
    result.readid = workItem.read.id;
	result.kd.setID(result.readid);
	
	std::string readSeq = workItem.read.seq.toString();	
	SeedVector seedVec, pieceVec;
	searchSeedsWithHybridKmers(readSeq, seedVec, result);
	//Push the first seed into pieceVec, which will be popped later as source seed
	if	(
		seedVec.size() >= 2
		&& !m_params.OnlySeed
		)
	{
		result.correctedLen += seedVec[0].seedStr.length();		
		pieceVec.push_back(seedVec[0]);		
	}
	else
		//Give up reads with less than 2 seeds
		return result;
    //Reserve sufficient str length for fast append
    pieceVec.back().seedStr.reserve(readSeq.length());
    initCorrect(readSeq, seedVec, pieceVec, result);
	
	result.merge = true;
	result.totalReadsLen = readSeq.length();
	for(SeedVector::iterator iter = pieceVec.begin(); iter != pieceVec.end(); iter++)
		result.correctedStrs.push_back((*iter).seedStr);
	
	return result;
    
}
void PacBioSelfCorrectionProcess::initCorrect(std::string& readSeq, const SeedVector& seedVec, SeedVector& pieceVec, PacBioSelfCorrectionResult& result)
{
	
	for(SeedVector::const_iterator iterTarget = seedVec.begin() + 1; iterTarget != seedVec.end(); iterTarget++)
	{
		int isFMExtensionSuccess = 0, firstFMExtensionType = 0;
		SeedFeature& source = pieceVec.back();
		std::string mergedSeq;
		
		for(int next = 0; next < m_params.numOfNextTarget && (iterTarget + next) != seedVec.end() ; next++)
		{
			const SeedFeature& target = *(iterTarget + next);
			isFMExtensionSuccess = 
			correctByFMExtension(source, target, readSeq, mergedSeq, result);
			firstFMExtensionType = (next == 0 ? isFMExtensionSuccess : firstFMExtensionType);
			if(isFMExtensionSuccess > 0)
			{
				result.totalWalkNum++;
				source.append(mergedSeq, target);
				iterTarget += next;
				break;
			}
			
		}		
		
		if(isFMExtensionSuccess <= 0)
		{
			const SeedFeature& target = *iterTarget;
			switch(firstFMExtensionType)
			{
				case -1: 
					result.highErrorNum++;
					break;
				case -2:
					result.exceedDepthNum++;
					break;
				case -3:
					result.exceedLeaveNum++;
					break;
				case -4:
					result.highErrorNum++;
					break;
				default:
					std::cout << "Does it really happen?\n";
					exit(EXIT_FAILURE);
			}
			result.totalWalkNum++;
			bool isMSAlignmentSuccess = 
			correctByMSAlignment(source, target, readSeq, mergedSeq, result);
			if(isMSAlignmentSuccess)
				source.append(mergedSeq, target);
			else if(!m_params.isSplit)
			{
				mergedSeq = readSeq.substr(source.seedEndPos + 1, target.seedEndPos - source.seedEndPos);
				source.append(mergedSeq, target);
				result.correctedLen += target.seedStr.length();
			}
			else
			{
				pieceVec.push_back(target);
				result.correctedLen += target.seedStr.length();
			}
		}
	
	}
	
	/*
	for(size_t targetSeed = 1 ; targetSeed < seedVec.size() ; targetSeed++)
	{	
		int FMWalkReturnType = 0;
		SeedFeature source = pieceVec.back();
		SeedFeature target =  seedVec[targetSeed];
		size_t extendKmerSize = std::min(source.endBestKmerSize, seedVec[targetSeed].startBestKmerSize) - 2;
			
        
		for(int nextTargetSeed = 0 ; nextTargetSeed < (int)m_params.numOfNextTarget && targetSeed + nextTargetSeed < seedVec.size() ; nextTargetSeed++)
		{
			target =  seedVec[targetSeed + nextTargetSeed];
			if((source.isRepeat || target.isRepeat) )
			{
				extendKmerSize = std::min(source.seedLength, target.seedLength);
				if(extendKmerSize > m_params.kmerLength+2) 
						extendKmerSize = m_params.kmerLength+2;
			}

			//int dis_between_src_target = target.seedStartPos - seedVec[targetSeed-1].seedStartPos - seedVec[targetSeed-1].seedStr.length();
			int dis_between_src_target = target.seedStartPos - seedVec[targetSeed-1].seedEndPos - 1;
			std::string mergedseq;
			FMWalkReturnType = extendBetweenSeeds(source, target, readSeq, mergedseq, extendKmerSize, dis_between_src_target, result);
			if(FMWalkReturnType > 0)
			{
				pieceVec.back().append(mergedseq);
				pieceVec.back().endBestKmerSize = target.endBestKmerSize;
				pieceVec.back().isRepeat = target.isRepeat;
                pieceVec.back().maxFixedMerFreqs = target.maxFixedMerFreqs;
				result.correctedLen += mergedseq.length();
                if(FMWalkReturnType ==1)
                    result.FMNum++;
                else
                    result.DPNum++;
				result.seedDis += dis_between_src_target;
				
				targetSeed = targetSeed + nextTargetSeed;
                result.totalWalkNum++;
				break;
			}
		}
		
		// All targets failure: 
		// 0: seed inter-distance too large
		// -1: kmer extension failed at later stage close to target
		// -4: kmer extension failed at early stage close to source
		// -2: exceed depth
		// -3: exceed leaves
		if(FMWalkReturnType <= 0)
		{
			//result.seedDis += seedVec[targetSeed].seedStartPos - seedVec[targetSeed-1].seedStartPos - seedVec[targetSeed-1].seedStr.length();
			result.seedDis += seedVec[targetSeed].seedStartPos - seedVec[targetSeed-1].seedEndPos - 1;
			result.correctedLen += seedVec[targetSeed].seedStr.length();

			if(!m_params.isSplit)
			{
				//size_t startPos = seedVec[targetSeed-1].seedStartPos + seedVec[targetSeed-1].seedStr.length();
				size_t startPos = seedVec[targetSeed-1].seedEndPos + 1;
				//size_t extendedLen = seedVec[targetSeed].seedStartPos + seedVec[targetSeed].seedStr.length() - startPos;
				size_t extendedLen = seedVec[targetSeed].seedEndPos - seedVec[targetSeed - 1].seedEndPos;
				pieceVec.back().append(readSeq.substr(startPos,extendedLen));
				pieceVec.back().endBestKmerSize = seedVec[targetSeed].endBestKmerSize;
				pieceVec.back().isRepeat = seedVec[targetSeed].isRepeat;
                pieceVec.back().maxFixedMerFreqs = target.maxFixedMerFreqs;
			}
			else
			{
				pieceVec.push_back(seedVec[targetSeed]);
			}
			
			if(FMWalkReturnType == -1 || FMWalkReturnType == -4)
				result.highErrorNum++;
			else if(FMWalkReturnType == -2)
				result.exceedDepthNum++;
			else if(FMWalkReturnType == -3)
				result.exceedLeaveNum++;
            if(FMWalkReturnType != 0)
            result.totalWalkNum++;
		}
	}
	*/
}

// Search seeds with fixed and dynamic kmer size. Noted by KuanWeiLee 20171027
void PacBioSelfCorrectionProcess::searchSeedsWithHybridKmers(const std::string& readSeq, SeedVector& seedVec, PacBioSelfCorrectionResult& result)
{
    Timer* seedTimer = new Timer("Seed Time",true);
    const unsigned int staticKmerSize = m_params.kmerLength;
	if(readSeq.length() < staticKmerSize) return;
    
	//Part 1 : Slide with a fixed kmer. Noted by KuanWeiLee
	std::vector<BistrandBWTInterval> FixedMerInterval;
	
    for(size_t i = 0 ; i <= readSeq.length() - staticKmerSize ; i++)
	{
        std::string kmer = readSeq.substr(i,staticKmerSize);
		BistrandBWTInterval object;
		object.fwdInterval = BWTAlgorithms::findInterval(m_params.indices.pRBWT, reverse(kmer));
		object.rvcInterval = BWTAlgorithms::findInterval(m_params.indices.pBWT, reverseComplement(kmer));
		FixedMerInterval.push_back(object);
		size_t currentKmerFreqs = object.getFreqs();
		result.kd.add(currentKmerFreqs);
    }
	result.kd.computeKDAttributes(1);
	//LOWCOV:1 UNIQUE:2 REPEAT:3
    KmerThresholdTable::TYPE mode = KmerThresholdTable::TYPE::UNIQUE;
	//KmerThresholdTable::TYPE mode = KmerThresholdTable::TYPE::REPEAT;
	
	float staticKmerThresholdValue, repeatKmerThresholdValue;
	staticKmerThresholdValue = KmerThresholdTable::get(staticKmerSize,mode);
	repeatKmerThresholdValue = mode == 3 ? staticKmerThresholdValue : staticKmerThresholdValue * 5;
	
	//Part 2 : Search seeds; slide through the read sequence with dynamic kmers. Noted by KuanWeiLee
	for(size_t startPos = 0; startPos <= readSeq.length() - staticKmerSize; startPos++)
	{
		bool isSeed = false, isRepeat = false;
		std::string kmer = readSeq.substr(startPos,staticKmerSize);
		unsigned int dynamicKmerSize = staticKmerSize;
		BWTInterval fwdInterval = FixedMerInterval[startPos].fwdInterval,
					rvcInterval = FixedMerInterval[startPos].rvcInterval;
		size_t fwdKmerFreqs, rvcKmerFreqs, currentKmerFreqs;
		size_t fixedMerFreqs, maxFixedMerFreqs = FixedMerInterval[startPos].getFreqs();
		float dynamicKmerThresholdValue;
		size_t seedStartPos = startPos, seedEndPos;
		float GCratio, freqsDiff;
		for(size_t movePos = startPos; movePos <= readSeq.length() - staticKmerSize; movePos++, isSeed = true)
		{	
			
			if(isSeed)
			{
				char s = readSeq[movePos + staticKmerSize - 1];//s:step
				char rcs = complement(s);
				kmer += s;
				dynamicKmerSize++;
				BWTAlgorithms::updateInterval(fwdInterval,s,m_params.indices.pRBWT);
				BWTAlgorithms::updateInterval(rvcInterval,rcs,m_params.indices.pBWT);
			}
			isRepeat = isRepeat || (maxFixedMerFreqs >= repeatKmerThresholdValue);
			dynamicKmerThresholdValue = KmerThresholdTable::get(dynamicKmerSize,mode);
			fwdKmerFreqs = fwdInterval.getFreqs();
			rvcKmerFreqs = rvcInterval.getFreqs();			
			currentKmerFreqs = fwdKmerFreqs + rvcKmerFreqs;
			fixedMerFreqs = FixedMerInterval[movePos].getFreqs();
			freqsDiff = (float)fixedMerFreqs / (float)maxFixedMerFreqs;
			//Gerneral seed extension strategy. Noted by KuanWeiLee
			if  (
				   dynamicKmerSize > m_params.kmerLengthUpperBound											//1.over length
				|| isLowComplexity(kmer,GCratio)															//2.low complexity
				|| fixedMerFreqs < staticKmerThresholdValue													//3.fixed kmer frequency
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
			//The order between general seed extentsion and kmer hitchhike strategy may make difference 
			//,which needs furher observation. Noted by KuanWeiLee 20171027
			//Kmer Hitchhike strategy.
			if	(isRepeat && freqsDiff < m_params.khhRatio)													//5.hitchhiking kmer (HIGH-->LOW)
			{
				kmer.erase(--dynamicKmerSize);
				seedEndPos = seedStartPos + dynamicKmerSize - 1;
				startPos = seedEndPos + 1;
				break;
			}
			if	(isRepeat && freqsDiff > m_params.r_khhRatio)												//6.hitchhiking kmer (LOW-->HIGH)
			{
				isSeed = false;
				startPos = movePos - 1;
				break;
			}
			
			maxFixedMerFreqs = MAX(maxFixedMerFreqs,fixedMerFreqs);
		}
		if(isSeed)
		{
			SeedFeature newSeed(seedStartPos, kmer, isRepeat, staticKmerSize, m_params.PBcoverage>>1, maxFixedMerFreqs);
			newSeed.estimateBestKmerSize(m_params.indices);
			seedVec.push_back(newSeed);
		}
	}
	seedVec = removeHitchhikingSeeds(seedVec, result);
	//[Debugseed] Output seeds & kd statistics for each reads.Noted by KuanWeiLee
    if(m_params.DebugSeed)
    {
        std::ostream* pSeedWriter = createWriter(m_params.directory + "seed/" + result.readid + ".seed");
		write(*pSeedWriter, seedVec);
		delete pSeedWriter;
		
		std::ostream* pKdWriter = createWriter(m_params.directory + "kmer-stat", std::ios_base::app);
		result.kd.write(*pKdWriter,KmerDistribution::TYPE::ATTRIBUTE);
		delete pKdWriter;
	}
	
	result.totalSeedNum = seedVec.size();
	result.Timer_Seed = seedTimer->getElapsedWallTime(); 
    delete seedTimer;
}

PacBioSelfCorrectionProcess::SeedVector PacBioSelfCorrectionProcess::removeHitchhikingSeeds(SeedVector initSeedVec, PacBioSelfCorrectionResult& result)
{
	for(SeedVector::iterator iterQuery = initSeedVec.begin(); iterQuery != initSeedVec.end(); iterQuery++)
	{
		SeedFeature& query = *iterQuery;
		SeedVector::iterator iterTarget = iterQuery + 1;
		//if(query.isHitchhiked) continue;
		while(iterTarget != initSeedVec.end() && (*iterTarget).seedStartPos - query.seedEndPos <= m_params.repaetDistance)
		{
			SeedFeature& target = *iterTarget;
			iterTarget++;
			//if(target.isHitchhiked) continue;
			float freqsDiff = (float)target.maxFixedMerFreqs / (float)query.maxFixedMerFreqs;
			target.isHitchhiked = target.isHitchhiked || (query.isRepeat && freqsDiff < m_params.shhRatio); //HIGH --> LOW
			query.isHitchhiked = query.isHitchhiked || (target.isRepeat && freqsDiff > m_params.r_shhRatio);//LOW  --> HIGH
		}
	}
	unsigned int index = 0;
	SeedVector finalSeedVec, outcastSeedVec;
	finalSeedVec.reserve(initSeedVec.size());
	outcastSeedVec.reserve(initSeedVec.size()>>3);
	for(SeedVector::iterator iter = initSeedVec.begin(); iter != initSeedVec.end(); iter++, index++)
	{
		if((*iter).isHitchhiked)
			outcastSeedVec.push_back(*iter);
		else
			finalSeedVec.push_back(*iter);
	}
	if(m_params.DebugSeed)
	{
		std::ostream* pOutcastSeedWriter = createWriter(m_params.directory + "seed/shh/" + result.readid + ".seed");
		write(*pOutcastSeedWriter, outcastSeedVec);
		delete pOutcastSeedWriter;
	}
	return finalSeedVec;
}

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

void PacBioSelfCorrectionProcess::write(std::ostream& outfile, const SeedVector& seedVec) const
{
	for(SeedVector::const_iterator iter = seedVec.begin(); iter != seedVec.end(); iter++)
		outfile << (*iter).seedStr << "\t" << (*iter).maxFixedMerFreqs << "\t" 
		<< (*iter).seedStartPos << "\t" << ((*iter).isRepeat ? "Yes" : "No") << "\n";
}

int PacBioSelfCorrectionProcess::correctByFMExtension
(const SeedFeature& source, const SeedFeature& target, const std::string& in, std::string& out, PacBioSelfCorrectionResult& result)
{
	Timer* FMTimer = new Timer("FM Time",true);
	int interval = target.seedStartPos - source.seedEndPos - 1;
	size_t extendKmerSize = MIN(source.endBestKmerSize, target.startBestKmerSize) - 2;
	if(source.isRepeat || target.isRepeat)
	{
		extendKmerSize = MIN(source.seedLength, target.seedLength);
		extendKmerSize = MIN(extendKmerSize, m_params.kmerLength + 2);
	}
	std::string src, trg, path;
	src = source.seedStr.substr(source.seedLength - extendKmerSize);
	trg = target.seedStr;
	path = in.substr(source.seedEndPos + 1, interval);
	int min_SA_threshold = 3, isFMExtensionSuccess = 0;
	min_SA_threshold = m_params.PBcoverage > 60 ? ((m_params.PBcoverage / 60) * 3) : min_SA_threshold;
	//float freqsDiff = (float)target.maxFixedMerFreqs / (float)source.maxFixedMerFreqs;
	bool isFromRtoU = false;
	//isFromRtoU = source.isRepeat && freqsDiff < m_params.shhRatio;//HIGH --> LOW
	isFromRtoU = source.isRepeat && source.maxFixedMerFreqs > target.maxFixedMerFreqs;
	if(isFromRtoU)
	{
		std::swap(src,trg);
		src = reverseComplement(src);
		trg = reverseComplement(trg);
		path = reverseComplement(path);
	}
	
	FMWalkResult2 fmwalkresult;
	LongReadSelfCorrectByOverlap OverlapTree
	(src, path, trg, interval, extendKmerSize, extendKmerSize + 2, m_params.FM_params, min_SA_threshold);
	isFMExtensionSuccess = OverlapTree.extendOverlap(fmwalkresult);
	result.Timer_FM += FMTimer->getElapsedWallTime();
	delete FMTimer;
	
	if(isFMExtensionSuccess < 0) return isFMExtensionSuccess;
	if(isFromRtoU)
	{
		fmwalkresult.mergedSeq = reverseComplement(fmwalkresult.mergedSeq);
		fmwalkresult.mergedSeq += reverseComplement(src).substr(extendKmerSize);
	}
	out = fmwalkresult.mergedSeq;
	out.erase(0,extendKmerSize);
	result.correctedLen += out.length();
	result.seedDis += interval;
	result.FMNum++;
	return isFMExtensionSuccess;
}

bool PacBioSelfCorrectionProcess::correctByMSAlignment
(const SeedFeature& source, const SeedFeature& target, const std::string& in, std::string& out, PacBioSelfCorrectionResult& result)
{
	if(m_params.NoDp) return false;
	Timer* DPTimer = new Timer("DP Time", true);
	int interval = target.seedStartPos - source.seedEndPos - 1;
	size_t extendKmerSize = MIN(source.endBestKmerSize, target.startBestKmerSize) - 2;
	if(source.isRepeat || target.isRepeat)
	{
		extendKmerSize = MIN(source.seedLength, target.seedLength);
		extendKmerSize = MIN(extendKmerSize, m_params.kmerLength + 2);
	}
	std::string src, trg, path;
	src = source.seedStr.substr(source.seedLength - extendKmerSize);
	trg = target.seedStr;
	path = in.substr(source.seedEndPos + 1, interval);
	path = src + path + trg;
	double identity = 0.65;
	size_t totalMaxFixedMerFreqs = source.maxFixedMerFreqs + target.maxFixedMerFreqs, min_call_coverage = 15;
	identity += (totalMaxFixedMerFreqs > 50  ? 0.5 : 0);
	identity += (totalMaxFixedMerFreqs > 100 ? 0.5 : 0);
	min_call_coverage = totalMaxFixedMerFreqs > 50 ? totalMaxFixedMerFreqs * 0.4 : min_call_coverage;
	
	
	MultipleAlignment maquery = 
	LongReadOverlap::buildMultipleAlignment
	(path, extendKmerSize, extendKmerSize, path.length()/10, identity, m_params.PBcoverage, m_params.indices);
	result.Timer_DP += DPTimer->getElapsedWallTime();
	delete DPTimer;
	
	
	if(maquery.getNumRows() <= 3) return false;
	out = maquery.calculateBaseConsensus(min_call_coverage, -1);
	out.erase(0,extendKmerSize);
	result.correctedLen += out.length();
	result.seedDis += interval;
	result.DPNum++;
	return true;
}

// Perform FMindex extension between source and target seeds
// Return FMWalkReturnType
/*
int PacBioSelfCorrectionProcess::extendBetweenSeeds
(const SeedFeature& source, const SeedFeature& target, std::string& rawSeq, std::string& mergedseq, size_t extendKmerSize, size_t dis_between_src_target, PacBioSelfCorrectionResult& result)
{
   // size_t srcKmerSize = std::max(source.endBestKmerSize, extendKmerSize);
   std::string srcStr, trgStr;
   srcStr = source.seedStr.substr(source.seedStr.length()-extendKmerSize);
   trgStr = target.seedStr;
   std::string strbetweensrctarget = rawSeq.substr(target.seedStartPos-dis_between_src_target,dis_between_src_target);

   int min_SA_threshold =3;
    //v1
    if (m_params.PBcoverage > 60) min_SA_threshold = int(m_params.PBcoverage/60)*3;

	
    
    FMWalkResult2 fmwalkresult;
    Timer* FMTimer = new Timer("FM Time",true);
    int FMWalkReturnType =0;
	if(source.isRepeat && source.maxFixedMerFreqs > target.maxFixedMerFreqs  )
	{
		LongReadSelfCorrectByOverlap OverlapTree(reverseComplement(trgStr),reverseComplement(strbetweensrctarget),reverseComplement(srcStr),dis_between_src_target,extendKmerSize,extendKmerSize+2,m_params.FM_params,min_SA_threshold);
		FMWalkReturnType = OverlapTree.extendOverlap(fmwalkresult);
		fmwalkresult.mergedSeq = reverseComplement(fmwalkresult.mergedSeq) + trgStr.substr(extendKmerSize);
  
	}
    else
    {
		LongReadSelfCorrectByOverlap OverlapTree(srcStr,strbetweensrctarget,trgStr,dis_between_src_target,extendKmerSize,extendKmerSize+2,m_params.FM_params,min_SA_threshold);
        FMWalkReturnType = 	OverlapTree.extendOverlap(fmwalkresult);  
      
    }
    
       if(FMWalkReturnType > 0)//extend success by fm extend
        {
            
          mergedseq = fmwalkresult.mergedSeq;
          //mergedseq = mergedseq.substr(extendKmerSize);
          mergedseq.erase(0,extendKmerSize);
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
		std::string rawSubseq = srcStr +  strbetweensrctarget + trgStr;
        MultipleAlignment maquery = LongReadOverlap::buildMultipleAlignment(rawSubseq,
                                                        extendKmerSize, //m_params.PBKmerLength, 
                                                        extendKmerSize, //m_params.PBKmerLength,
                                                        rawSubseq.length()/10, 
                                                        identity,	// alignment identity < 0.7 are often false positive repeats
                                                        m_params.PBcoverage,
                                                        m_params.indices);
        
        std::string consensus = maquery.calculateBaseConsensus(min_call_coverage, -1);
        if(maquery.getNumRows() <= 3)
                return FMWalkReturnType;
        mergedseq = consensus;
        if(!mergedseq.empty())
			mergedseq.erase(0,extendKmerSize);
        
        result.Timer_DP += DPTimer->getElapsedWallTime() ;
        delete DPTimer;
        
        FMWalkReturnType = 2;
    
    }
    
     
	return FMWalkReturnType;
}
*/
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
m_highErrorNum(0),
m_exceedDepthNum(0),
m_exceedLeaveNum(0),
m_FMNum(0),
m_DPNum(0),
m_OutcastNum(0),
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
	m_OutcastNum = m_totalWalkNum - m_FMNum - m_DPNum;
	if(m_totalWalkNum>0 && m_totalReadsLen>0)
	{
		std::cout << "\n"
		<< "totalReadsLen: " << m_totalReadsLen << "\n"
		<< "correctedLen: " << m_correctedLen << ", ratio: " << (float)(m_correctedLen)/m_totalReadsLen << "\n"
		<< "totalSeedNum: " << m_totalSeedNum << "\n"
		<< "totalWalkNum: " << m_totalWalkNum << "\n"
		<< "FMNum: " << m_FMNum << ", ratio: " << (float)(m_FMNum * 100)/m_totalWalkNum << "%\n"
		<< "DPNum: " << m_DPNum << ", ratio: " << (float)(m_DPNum * 100)/m_totalWalkNum << "%\n"
        << "OutcastNum: " << m_OutcastNum << ", ratio: " << (float)(m_OutcastNum * 100)/m_totalWalkNum << "%\n"
		<< "highErrorNum: " << m_highErrorNum << ", ratio: " << (float)(m_highErrorNum * 100)/(m_DPNum + m_OutcastNum) << "%\n"
		<< "exceedDepthNum: " << m_exceedDepthNum << ", ratio: " << (float)(m_exceedDepthNum * 100)/(m_DPNum + m_OutcastNum) << "%\n"
		<< "exceedLeaveNum: " << m_exceedLeaveNum << ", ratio: " << (float)(m_exceedLeaveNum * 100)/(m_DPNum + m_OutcastNum) << "%\n"
		<< "disBetweenSeeds: " << m_seedDis/m_totalWalkNum << "\n"
        << "Time of searching Seeds: " << m_Timer_Seed  << "\n"
        << "Time of searching FM: " << m_Timer_FM  << "\n"
        << "Time of searching DP: " << m_Timer_DP  << "\n";
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
		m_highErrorNum += result.highErrorNum;
		m_exceedDepthNum += result.exceedDepthNum;
		m_exceedLeaveNum += result.exceedLeaveNum;
		m_FMNum += result.FMNum;
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
		
		//for(size_t i = 0 ; i < result.correctedStrs.size() ; i++)
		for(std::vector<DNAString>::const_iterator iter = result.correctedStrs.begin(); iter != result.correctedStrs.end(); iter++)
		{
			size_t i = iter - result.correctedStrs.begin();
			SeqItem mergeRecord2;
			std::stringstream ss2;
			//ss2 << item.read.id << "_" << i << "_" << result.correctedStrs[i].toString().length();
			ss2 << item.read.id << "_" << i << (*iter).toString().length();
			mergeRecord2.id = ss2.str();
			//mergeRecord2.seq = result.correctedStrs[i];
			mergeRecord2.seq = *iter;
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