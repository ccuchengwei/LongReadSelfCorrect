//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioSelfCorrectionProcess - Self-correction using FM-index walk for PacBio reads
//

#ifndef PACBIOSELFCORRECTIONPROCESS_H
#define PACBIOSELFCORRECTIONPROCESS_H

#include "HashMap.h"
#include "Util.h"
#include "SequenceProcessFramework.h"
#include "SequenceWorkItem.h"
#include "Metrics.h"
#include "BWTIndexSet.h"
#include "SampledSuffixArray.h"
#include "BWTAlgorithms.h"
#include "KmerDistribution.h"
#include "LongReadCorrectByOverlap.h"
#include "SeedFeature.h"

// Parameter object for the error corrector
struct PacBioSelfCorrectionParameters
{
	//PacBioSelfCorrectionAlgorithm algorithm;
	BWTIndexSet indices;
	FMextendParameters FM_params;

	int numKmerRounds;
	unsigned int kmerLength;
	unsigned int kmerLengthUpperBound;
	//unsigned int repaetDistance = 50;//70-->100. Need further check. Noted by KuanWeiLee 20171125
	unsigned int repaetDistance = 100;
	//float shhRatio = 0.4f;
	float shhRatio = 0.6f;
	float r_shhRatio = 1.f/shhRatio;
	//float khhRatio = 0.4f;
	float khhRatio = 0.6f;
	float r_khhRatio = 1.f/khhRatio;
	//The hhRatio for seed or kmer hitchhike may be different; 
	//in repeat mode, one for seed seems lower but the other for kmer still remains unknown.Noted by KuanWeiLee 20171213
	
	// tree search parameters
	int maxLeaves;
	int minOverlap;
	int maxOverlap;
    int idmerLength;
    

	// PACBIO
	int minKmerLength;
	int FMWKmerThreshold;
	int seedKmerThreshold;
	int numOfNextTarget;
	int collectedSeeds;
    double ErrorRate;
	bool isSplit;
	bool isFirst;
	size_t maxSeedInterval;
	int PBcoverage;
    
    bool DebugExtend;
    bool DebugSeed;
	bool OnlySeed;
	bool NoDp;
	std::string directory;
};




struct PacBioSelfCorrectionResult
{
	PacBioSelfCorrectionResult()
	: merge(false),
	totalReadsLen(0),
	correctedLen(0),
	totalSeedNum(0),
	totalWalkNum(0),
	highErrorNum(0),
	exceedDepthNum(0),
	exceedLeaveNum(0),
	FMNum(0),
    DPNum(0),
	seedDis(0),
    Timer_Seed(0),
    Timer_FM(0),
    Timer_DP(0) {}

	std::string readid;
	KmerDistribution kd;
	
	bool merge;
	
	size_t kmerLength;

	// PacBio reads correction by Ya, v20151001.
	std::vector<DNAString> correctedStrs;
	int64_t totalReadsLen;
	int64_t correctedLen;
	int64_t totalSeedNum;
	int64_t totalWalkNum;
	int64_t highErrorNum;
	int64_t exceedDepthNum;
	int64_t exceedLeaveNum;
	int64_t FMNum;
    int64_t DPNum;
	int64_t seedDis;
    double Timer_Seed;
    double Timer_FM;
    double Timer_DP;
};

//
class PacBioSelfCorrectionProcess
{
public:

	PacBioSelfCorrectionProcess(const PacBioSelfCorrectionParameters params):m_params(params){};
	~PacBioSelfCorrectionProcess(){};
	PacBioSelfCorrectionResult process(const SequenceWorkItem& workItem);

private:
	typedef std::vector<SeedFeature> SeedVector;
	PacBioSelfCorrectionParameters m_params;
	
	//search seeds
    void searchSeedsWithHybridKmers(const std::string& readSeq, SeedVector& seedVec, PacBioSelfCorrectionResult &result);
	SeedVector removeHitchhikingSeeds(SeedVector initSeedVec, PacBioSelfCorrectionResult& result);
	bool isLowComplexity (const std::string& seq, float & GCratio, float threshold = 0.7f);
	void write(std::ostream& outfile, const SeedVector& seedVec) const;
	//correct sequence
	void initCorrect(std::string& readSeq, const SeedVector& seeds, SeedVector& pacbioCorrectedStrs, PacBioSelfCorrectionResult& result);
	int correctByFMExtension
	(const SeedFeature& source, const SeedFeature& target, const std::string& in, std::string& out, PacBioSelfCorrectionResult& result);
	bool correctByMSAlignment
	(const SeedFeature& source, const SeedFeature& target, const std::string& in, std::string& out, PacBioSelfCorrectionResult& result);
	
};

//postprocess
class PacBioSelfCorrectionPostProcess
{
public:
	PacBioSelfCorrectionPostProcess(std::string correctFile, std::string discardFile, const PacBioSelfCorrectionParameters params);
	~PacBioSelfCorrectionPostProcess();
	void process(const SequenceWorkItem& item, const PacBioSelfCorrectionResult& result);
	
private:

	std::ostream* m_pCorrectWriter;
	std::ostream* m_pDiscardWriter;
	std::ostream* m_pKdWriter;
	PacBioSelfCorrectionParameters m_params;
	KmerDistribution m_kd;
	
	int64_t m_totalReadsLen;
	int64_t m_correctedLen;
	int64_t m_totalSeedNum;
	int64_t m_totalWalkNum;
	int64_t m_highErrorNum;
	int64_t m_exceedDepthNum;
	int64_t m_exceedLeaveNum;
	int64_t m_FMNum;
    int64_t m_DPNum;
	int64_t m_OutcastNum;
	int64_t m_seedDis;
    double m_Timer_Seed;
    double m_Timer_FM;
    double m_Timer_DP;

};


#endif
