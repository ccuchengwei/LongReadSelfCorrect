///-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioSelfCorrectionProcess.cpp - Correction of PacBio reads using FM-index
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include "Util.h"
#include "PacBioSelfCorrection.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "OverlapCommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "ASQG.h"
#include "gzstream.h"
#include "SequenceProcessFramework.h"
#include "PacBioSelfCorrectionProcess.h"
#include "CorrectionThresholds.h"
#include "BWTIntervalCache.h"
#include "KmerThresholdTable.h"

#define FORMULA( x,y,z ) ( (x) ? (0.05776992234f * y - 0.4583043394f * z + 10.19159685f) : (0.0710704607f * y - 0.5445663957f * z + 12.26253388f) )
//
// Getopt
//
#define SUBPROGRAM "PacBioSelfCorrection"
static const char *CORRECT_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"Written by Yao-Ting Huang & Ping-Yeh Chen.\n"
"\n"
"Copyright 2015 National Chung Cheng University\n";

static const char *CORRECT_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Correct PacBio reads via FM-index walk\n"
"\n"
"      --help                           Display this help and exit\n"
"      -v, --verbose                    Display verbose output\n"
"      -p, --prefix=PREFIX              Use PREFIX for the names of the index files\n"
"      -o, --directory=PATH             Put results in the directory\n"
"      -t, --threads=NUM                Use NUM threads for the computation (default: 1)\n"
"\nPacBio correction parameters:\n"
"      -k, --kmer-size=N                The length of the kmer to use. (default: 19 (PacBioS).)\n"
"      -s, --min-kmer-size=N            The minimum length of the kmer to use. (default: 13.)\n"
"      -x, --kmer-threshold=N           Attempt to correct kmers that are seen less than N times. (default: 3)\n"
"      -e, --error-rate=N               The error rate of PacBio reads.(default:0.15)\n"
"      -i, --idmer-length=N             The length of the kmer to identify similar reads.(default: 9)\n"
"      -L, --max-leaves=N               Number of maximum leaves in the search tree. (default: 32)\n"
"      -C, --PBcoverage=N               Coverage of PacBio reads(default: 90)\n"
"      --debugseed                      Output seeds file for each reads (default: false)\n"
"      --debugextend                    Show extension information (default: false)\n"
"      --onlyseed                       Only search seeds file for each reads (default: false)\n"
"      --nodp                           Don't use dp (default: false)\n"
"      --split                          Split the uncorrected reads (default: false)\n"

"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
	static unsigned int verbose;
	static int numThreads = 1;
	static std::string prefix;
	static std::string readsFile;
	static std::string correctFile;
	static std::string discardFile;
	static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
	static int kmerLength = 19;
	static int kmerThreshold = 3;
	static int maxLeaves=32;
    static int idmerLength = 9;
    static double ErrorRate=0.15;	
	static int minKmerLength = 13;	
	static int numOfNextTarget = 1;
	static int collect = 5;
	static int kmerLengthUpperBound = 50;
	
	static bool split = false;
	static bool isFirst = false;
	size_t maxSeedInterval = 500;
    static size_t PBcoverage = 90;// PB seed searh depth
    static bool DebugExtend = false;
    static bool DebugSeed = false;
	static bool OnlySeed = false;
	static bool NoDp = false;
	static std::string directory;
}

static const char* shortopts = "p:t:o:k:x:L:s:d:c:C:e:i:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_DISCARD, OPT_SPLIT, OPT_FIRST,OPT_DEBUGEXTEND,OPT_DEBUGSEED,OPT_ONLYSEED,OPT_NODP };

static const struct option longopts[] = {
	{ "threads",       required_argument, NULL, 't' },
	{ "directory",     required_argument, NULL, 'o' },
	{ "prefix",        required_argument, NULL, 'p' },
	{ "kmer-size",     required_argument, NULL, 'k' },
	{ "kmer-threshold",required_argument, NULL, 'x' },
	{ "max-leaves",    required_argument, NULL, 'L' },
	{ "min-kmer-size", required_argument, NULL, 's' },
    { "error-rate",    required_argument, NULL, 'e' },
    { "idmer-length",  required_argument, NULL, 'i' },
	{ "downward",      required_argument, NULL, 'd' },
	{ "collect",       required_argument, NULL, 'c' },
    { "PBcoverage",    required_argument, NULL, 'C' },
	{ "verbose",       no_argument,       NULL, 'v' },
	{ "split",         no_argument,       NULL, OPT_SPLIT },
	{ "first",         no_argument,       NULL, OPT_FIRST },
    { "debugextend",   no_argument,       NULL, OPT_DEBUGEXTEND },
    { "debugseed",     no_argument,       NULL, OPT_DEBUGSEED },
	{ "onlyseed",      no_argument,       NULL, OPT_ONLYSEED },
	{ "nodp",          no_argument,       NULL, OPT_NODP },
	{ "discard",       no_argument,       NULL, OPT_DISCARD },
	{ "help",          no_argument,       NULL, OPT_HELP },
	{ "version",       no_argument,       NULL, OPT_VERSION },
	{ NULL, 0, NULL, 0 }
};

//
// Main
//
int PacBioSelfCorrectionMain(int argc, char** argv)
{
	parsePacBioSelfCorrectionOptions(argc, argv);

	// Set the error correction parameters
	PacBioSelfCorrectionParameters ecParams;
	BWT *pBWT, *pRBWT;
	SampledSuffixArray* pSSA;

	// Load indices
	#pragma omp parallel
	{
		#pragma omp single nowait
		{	//Initialization of large BWT takes some time, pass the disk to next job
			std::cout << "Loading BWT: " << opt::prefix + BWT_EXT << "\n";
			pBWT = new BWT(opt::prefix + BWT_EXT, opt::sampleRate);
		}
		#pragma omp single nowait
		{
			std::cout << "Loading RBWT: " << opt::prefix + RBWT_EXT << "\n";
			pRBWT = new BWT(opt::prefix + RBWT_EXT, opt::sampleRate);
		}
		#pragma omp single nowait
		{
			std::cout << "Loading Sampled Suffix Array: " << opt::prefix + SAI_EXT << "\n";
			pSSA = new SampledSuffixArray(opt::prefix + SAI_EXT, SSA_FT_SAI);
		}
	}
	BWTIndexSet indexSet;
	indexSet.pBWT = pBWT;
	indexSet.pRBWT = pRBWT;
	indexSet.pSSA = pSSA;
	ecParams.indices = indexSet;
	
	KmerThresholdTable::m_startLen = opt::kmerLength;
	KmerThresholdTable::m_endLen = opt::kmerLengthUpperBound;
	KmerThresholdTable::m_coverage = opt::PBcoverage;
	KmerThresholdTable::m_lowcov = new float[opt::kmerLengthUpperBound + 1]{};
	KmerThresholdTable::m_unique = new float[opt::kmerLengthUpperBound + 1]{};
	KmerThresholdTable::m_repeat = new float[opt::kmerLengthUpperBound + 1]{};
	KmerThresholdTable::pTableWriter = createWriter(opt::directory + "threshold-table");	
	KmerThresholdTable::compute();
	KmerThresholdTable::write();

	
	// Open outfiles and start a timer
	Timer* pTimer = new Timer(PROGRAM_IDENT);

	ecParams.kmerLength = opt::kmerLength;
	ecParams.kmerLengthUpperBound = opt::kmerLengthUpperBound;
	ecParams.maxLeaves = opt::maxLeaves;
	ecParams.minKmerLength = opt::minKmerLength;
    ecParams.idmerLength = opt::idmerLength;
    ecParams.ErrorRate = opt::ErrorRate;
	ecParams.FMWKmerThreshold = opt::kmerThreshold;
	ecParams.numOfNextTarget = opt::numOfNextTarget;
	ecParams.collectedSeeds = opt::collect;
    ecParams.PBcoverage = opt::PBcoverage;
	ecParams.isSplit = opt::split;
	ecParams.isFirst = opt::isFirst;
    ecParams.DebugExtend = opt::DebugExtend;
    ecParams.DebugSeed = opt::DebugSeed;
	ecParams.OnlySeed = opt::OnlySeed;
	ecParams.NoDp = opt::NoDp;
	ecParams.maxSeedInterval = opt::maxSeedInterval;
	ecParams.directory = opt::directory;
	

	std::cout << "\nCorrecting PacBio reads for " << opt::readsFile << " using--\n"
	<< "number of threads:\t" << opt::numThreads << "\n"
	<< "PB reads coverage:\t" << ecParams.PBcoverage << "\n"
	<< "large kmer size:\t" << ecParams.kmerLength << "\n" 
	<< "small kmer size:\t" << ecParams.minKmerLength << "\n"
	<< "small kmer freq. cutoff:\t" << ecParams.FMWKmerThreshold << "\n"
	<< "max leaves:\t" << ecParams.maxLeaves  << "\n"
	<< "max depth:\t1.2~0.8* (length between two seeds +- 20)" << "\n"
	<< "num of next Targets:\t" << ecParams.numOfNextTarget << "\n";

	
	// Setup post-processor
	PacBioSelfCorrectionPostProcess* pPostProcessor = new PacBioSelfCorrectionPostProcess(opt::correctFile, opt::discardFile, ecParams);

	if(opt::numThreads <= 1)
	{
		// Serial mode
		PacBioSelfCorrectionProcess* pProcessor = new PacBioSelfCorrectionProcess(ecParams);

		SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
		PacBioSelfCorrectionResult,
		PacBioSelfCorrectionProcess,
		PacBioSelfCorrectionPostProcess>(opt::readsFile, pProcessor, pPostProcessor);
		delete pProcessor;
	}
	else
	{
		// Parallel mode
		std::vector<PacBioSelfCorrectionProcess*> pProcessorVector;
		for(int i = 0; i < opt::numThreads; ++i)
		{
			PacBioSelfCorrectionProcess* pProcessor = new PacBioSelfCorrectionProcess(ecParams);
			pProcessorVector.push_back(pProcessor);
		}

		SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
		PacBioSelfCorrectionResult,
		PacBioSelfCorrectionProcess,
		PacBioSelfCorrectionPostProcess>(opt::readsFile, pProcessorVector, pPostProcessor);

		// SequenceProcessFramework::processSequencesParallelOpenMP<SequenceWorkItem,
		// PacBioSelfCorrectionResult,
		// PacBioSelfCorrectionProcess,
		// PacBioSelfCorrectionPostProcess>(opt::readsFile, pProcessorVector, pPostProcessor);
		
		while(!pProcessorVector.empty())
		{
			delete pProcessorVector.back();
			pProcessorVector.pop_back();
		}
	}
	delete pTimer;
	delete pPostProcessor;
	delete pBWT;
	delete pRBWT;
	delete pSSA;
	
	KmerThresholdTable::release();
	return 0;
}


//
// Handle command line arguments
//
void parsePacBioSelfCorrectionOptions(int argc, char** argv)
{
	optind = 1;	//reset getopt
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
	{
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c)
		{
		case 'p': arg >> opt::prefix; break;
		case 'o': arg >> opt::directory; break;
		case 't': arg >> opt::numThreads; break;
		case 'k': arg >> opt::kmerLength; break;
		case 'x': arg >> opt::kmerThreshold; break;
		case 'L': arg >> opt::maxLeaves; break;
		case 's': arg >> opt::minKmerLength; break;
        case 'e': arg >> opt::ErrorRate; break;
        case 'i': arg >> opt::idmerLength; break;
		case 'd': arg >> opt::numOfNextTarget; break;
		case 'c': arg >> opt::collect; break;
        case 'C': arg >> opt::PBcoverage; break;
		case 'v': opt::verbose++; break;
		case '?': die = true; break;
		case OPT_SPLIT: opt::split = true; break;
		case OPT_FIRST: opt::isFirst = true; break;
        case OPT_DEBUGEXTEND: opt::DebugExtend = true; break;
        case OPT_DEBUGSEED: opt::DebugSeed = true; break;
		case OPT_ONLYSEED: opt::OnlySeed = true; break;
		case OPT_NODP: opt::NoDp = true; break;
		case OPT_HELP:
			std::cout << CORRECT_USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		case OPT_VERSION:
			std::cout << CORRECT_VERSION_MESSAGE;
			exit(EXIT_SUCCESS);
		}
	}

	if (argc - optind < 1)
	{
		std::cerr << SUBPROGRAM ": missing arguments\n";
		die = true;
	}
	else if (argc - optind > 1)
	{
		std::cerr << SUBPROGRAM ": too many arguments\n";
		die = true;
	}

	if(opt::numThreads <= 0)
	{
		std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
		die = true;
	}


	if(opt::kmerLength <= 0)
	{
		std::cerr << SUBPROGRAM ": invalid kmer length: " << opt::kmerLength << ", must be greater than zero\n";
		die = true;
	}

	if(opt::kmerThreshold <= 0)
	{
		std::cerr << SUBPROGRAM ": invalid kmer threshold: " << opt::kmerThreshold << ", must be greater than zero\n";
		die = true;
	}
	CorrectionThresholds::Instance().setBaseMinSupport(opt::kmerThreshold);

	if(opt::prefix.empty())
	{
		std::cerr << SUBPROGRAM << ": no prefix\n";
		die = true;
	}
	
	if(opt::directory.empty())
	{
		std::cerr << SUBPROGRAM << ": no directory\n";
		die = true;
	}
	else
	{		
		opt::directory = opt::directory + "/";
		std::string workingDir = opt::directory + (opt::DebugSeed ? "seed/stat/" : "");
		if( system(("mkdir -p " + workingDir).c_str()) != 0)
		{
			std::cerr << SUBPROGRAM << ": something wrong in directory: " << opt::directory << "\n";
			die = true;
		}
	}
	if(die)
	{
		std::cout << "\n" << CORRECT_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}
	
	opt::readsFile = argv[optind++];
	std::string out_prefix = stripFilename(opt::readsFile);
	opt::correctFile = opt::directory + out_prefix + ".correct.fa";
	opt::discardFile = opt::directory + out_prefix + ".discard.fa";
	
	std::string outfilename = opt::directory + "threshold-table";
	std::ofstream outfile(outfilename);
	for(int i=opt::kmerLength; i<=50; i++)
	{
		float kmerThresholdValueWithLowCoverage = FORMULA(true,opt::PBcoverage,i);
		float kmerThresholdValue = FORMULA(false,opt::PBcoverage,i);
		outfile << i << "\t";
		outfile << (kmerThresholdValue < 5 ? 5 : kmerThresholdValue) << "\t";
		outfile << (kmerThresholdValueWithLowCoverage < 5 ? 5 : kmerThresholdValueWithLowCoverage) << "\n";
	}
	outfile.close();
}
