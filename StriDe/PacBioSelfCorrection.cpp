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
"      -t, --threads=NUM                Use NUM threads for the computation (default: 1)\n"
"      -p, --prefix=PREFIX              Use PREFIX for the names of the index files\n"
"      -o, --directory=PATH             Put results in the directory\n"
"\nPacBio correction parameters:\n"
"      -c, --PBcoverage=N               Coverage of PacBio reads (default: 90)\n"
"      -e, --error-rate=N               The error rate of PacBio reads.(default:0.15)\n"
"      -k, --kmer-size=N                The start kmer length (default: 19 (PacBioS).)\n"
"      -d, --num-of-next-target         The number of next FMWalk target seed(default: 1)\n"
"      -l, --max-leaves=N               Number of maximum leaves in the search tree. (default: 32)\n"
"      -i, --idmer-length=N             The length of the kmer to identify similar reads.(default: 9)\n"
"      -s, --min-kmer-size=N            The minimum length of the kmer to use. (default: 13.)\n"
"      -n, --genome=(5/10/100)[m]       Genome size of the species (default: 10m)\n"
"      -m, --mode=(0/1/2)               Mode in seed-searching (default: 1)\n"
"      -v, --verbose                    Display verbose output\n"
"      --help                           Display this help and exit\n"
"      --version                        Display version and exit\n"
"      --manual                         Disable auto detection module (default: false)\n"
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
	static int numThreads = 1;
	static std::string prefix;
	static std::string directory;
	static std::string readsFile;
	static std::string correctFile;
	static std::string discardFile;
	static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
	
    static size_t PBcoverage = 90;// PB seed searh depth
    static double ErrorRate=0.15;
	static int startKmerLength = 19;
	static int numOfNextTarget = 1;
	static int maxLeaves=32;
    static int idmerLength = 9;
	static int minKmerLength = 13;
	
	static int genome = 10;
	static int mode = 1;
	static int verbose = 0;
	
	static bool Manual = false;
	static bool Split = false;
    static bool DebugExtend = false;
    static bool DebugSeed = false;
	static bool OnlySeed = false;
	static bool NoDp = false;
}

static const char* shortopts = "t:p:o:c:e:k:d:l:i:s:n:m:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_MANUAL, OPT_SPLIT, OPT_FIRST, OPT_DEBUGEXTEND, OPT_DEBUGSEED, OPT_ONLYSEED, OPT_NODP };

static const struct option longopts[] = {
	{ "threads",            required_argument, nullptr, 't' },
	{ "prefix",             required_argument, nullptr, 'p' },
	{ "directory",          required_argument, nullptr, 'o' },
    { "PBcoverage",         required_argument, nullptr, 'c' },
    { "error-rate",         required_argument, nullptr, 'e' },
	{ "kmer-size",          required_argument, nullptr, 'k' },
	{ "num-of-next-target", required_argument, nullptr, 'd' },
	{ "max-leaves",         required_argument, nullptr, 'l' },
    { "idmer-length",       required_argument, nullptr, 'i' },
	{ "min-kmer-size",      required_argument, nullptr, 's' },
    { "genome",             required_argument, nullptr, 'n' },
    { "mode",               required_argument, nullptr, 'm' },
	{ "verbose",            no_argument,       nullptr, 'v' },
	{ "help",               no_argument,       nullptr, OPT_HELP },
	{ "version",            no_argument,       nullptr, OPT_VERSION },
	{ "manual",             no_argument,       nullptr, OPT_MANUAL },
	{ "split",              no_argument,       nullptr, OPT_SPLIT },
    { "debugextend",        no_argument,       nullptr, OPT_DEBUGEXTEND },
    { "debugseed",          no_argument,       nullptr, OPT_DEBUGSEED },
	{ "onlyseed",           no_argument,       nullptr, OPT_ONLYSEED },
	{ "nodp",               no_argument,       nullptr, OPT_NODP },
	{ nullptr, 0, nullptr, 0 }
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
			std::cerr << "Loading BWT: " << opt::prefix + BWT_EXT << "\n";
			pBWT = new BWT(opt::prefix + BWT_EXT, opt::sampleRate);
		}
		#pragma omp single nowait
		{
			std::cerr << "Loading RBWT: " << opt::prefix + RBWT_EXT << "\n";
			pRBWT = new BWT(opt::prefix + RBWT_EXT, opt::sampleRate);
		}
		#pragma omp single nowait
		{
			std::cerr << "Loading Sampled Suffix Array: " << opt::prefix + SAI_EXT << "\n";
			pSSA = new SampledSuffixArray(opt::prefix + SAI_EXT, SSA_FT_SAI);
		}
	}
	BWTIndexSet indexSet;
	indexSet.pBWT = pBWT;
	indexSet.pRBWT = pRBWT;
	indexSet.pSSA = pSSA;
	ecParams.indices = indexSet;
	ecParams.directory = opt::directory;

    ecParams.PBcoverage = opt::PBcoverage;
    ecParams.ErrorRate = opt::ErrorRate;
	ecParams.startKmerLength = opt::startKmerLength;
	ecParams.numOfNextTarget = opt::numOfNextTarget;
	ecParams.maxLeaves = opt::maxLeaves;
    ecParams.idmerLength = opt::idmerLength;
	ecParams.minKmerLength = opt::minKmerLength;
	
	ecParams.Manual = opt::Manual;
	ecParams.Split = opt::Split;
    ecParams.DebugExtend = opt::DebugExtend;
    ecParams.DebugSeed = opt::DebugSeed;
	ecParams.OnlySeed = opt::OnlySeed;
	ecParams.NoDp = opt::NoDp;
	
	if(!opt::Manual)
	{
		switch(opt::genome)
		{
			case 5:
				ecParams.startKmerLength = 17;
				ecParams.kmerOffset[2] = 2;
				break;
			case 10:
				ecParams.startKmerLength = 19;
				ecParams.kmerOffset[2] = 4;
				break;
			case 100:
				ecParams.startKmerLength = 21;
				ecParams.kmerOffset[2] = 6;
				break;
		}
	}
	ecParams.mode = opt::mode;
	
	ecParams.kmerSet.insert(ecParams.startKmerLength);
	ecParams.kmerSet.insert(ecParams.scanKmerLength);
	for(auto& iter : ecParams.kmerOffset)
		ecParams.kmerSet.insert(ecParams.startKmerLength - iter);
	for(auto& iter : ecParams.overlapKmerLength)
		ecParams.kmerSet.insert(iter);
	
	FMextendParameters FM_params(
			indexSet,opt::idmerLength,
			opt::maxLeaves,
			opt::minKmerLength,
			opt::PBcoverage,
			opt::ErrorRate,
			false);
		//	opt::DebugExtend);
	ecParams.FM_params = FM_params;
	
	//Initialize KmerThresholdTable
	KmerThresholdTable::initialize(
			ecParams.startKmerLength,
			ecParams.kmerLengthUpperBound,
			ecParams.PBcoverage,
			ecParams.directory);
	KmerThresholdTable::compute();
	KmerThresholdTable::write();
	
	
	std::cerr << "\nCorrecting PacBio reads for " << opt::readsFile << " using--\n"
	<< "number of threads:\t" << opt::numThreads << "\n"
	<< "PB reads coverage:\t" << ecParams.PBcoverage << "\n"
	<< "num of next Targets:\t" << ecParams.numOfNextTarget << "\n"
	<< "large kmer size:\t" << ecParams.startKmerLength << "\n" 
	<< "small kmer size:\t" << ecParams.minKmerLength << "\n"
	<< "max leaves:\t" << ecParams.maxLeaves  << "\n"
	<< "max depth:\t1.2~0.8* (length between two seeds +- 20)" << "\n";

	
	// Start a timer
	Timer* pTimer = new Timer(PROGRAM_IDENT);
	
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
			pProcessorVector.push_back(new PacBioSelfCorrectionProcess(ecParams));

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
	delete pPostProcessor;
	delete pBWT;
	delete pRBWT;
	delete pSSA;
	
	KmerThresholdTable::release();
	delete pTimer;
	return 0;
}


//
// Handle command line arguments
//
void parsePacBioSelfCorrectionOptions(int argc, char** argv)
{
	optind = 1;	//reset getopt
	bool die = false;
	for (char c; (c = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;)
	{
		std::istringstream arg(optarg != nullptr ? optarg : "");
		switch (c)
		{
			case 't': arg >> opt::numThreads; break;
			case 'p': arg >> opt::prefix; break;
			case 'o': arg >> opt::directory; break;
			case 'c': arg >> opt::PBcoverage; break;
			case 'e': arg >> opt::ErrorRate; break;
			case 'k': arg >> opt::startKmerLength; break;
			case 'd': arg >> opt::numOfNextTarget; break;
			case 'l': arg >> opt::maxLeaves; break;
			case 'i': arg >> opt::idmerLength; break;
			case 's': arg >> opt::minKmerLength; break;
			case 'n': arg >> opt::genome; break;
			case 'm': arg >> opt::mode; break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				std::cerr << CORRECT_USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cerr << CORRECT_VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_MANUAL: opt::Manual = true; break;
			case OPT_SPLIT: opt::Split = true; break;
			case OPT_DEBUGEXTEND: opt::DebugExtend = true; break;
			case OPT_DEBUGSEED: opt::DebugSeed = true; break;
			case OPT_ONLYSEED: opt::OnlySeed = true; break;
			case OPT_NODP: opt::NoDp = true; break;
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
		opt::directory += "/";
		std::string workingDir = opt::directory + (opt::DebugSeed ? "seed/error/" : "");
		if( system(("mkdir -p " + workingDir).c_str()) != 0)
		{
			std::cerr << SUBPROGRAM << ": something wrong in directory: " << opt::directory << "\n";
			die = true;
		}
		workingDir = opt::directory + (opt::DebugSeed ? "extend/" : "");
		if( system(("mkdir -p " + workingDir).c_str()) != 0)
		{
			std::cerr << SUBPROGRAM << ": something wrong in directory: " << opt::directory << "\n";
			die = true;
		}
	}

	if(opt::PBcoverage <=0)
	{
		std::cerr << SUBPROGRAM ": invalid number of coverage: " << opt::PBcoverage << ", must be greater than zero\n";
		die = true;
	}
	
	if(opt::ErrorRate < 0 || opt::ErrorRate > 1)
	{
		std::cerr <<SUBPROGRAM ":invalid error rate: " << opt::ErrorRate << ", must be 0 ~ 1\n";
		die = true;
	}
	
	if(opt::startKmerLength <= 0)
	{
		std::cerr << SUBPROGRAM ": invalid start kmer length: " << opt::startKmerLength << ", must be greater than zero\n";
		die = true;
	}
	
	if(opt::numOfNextTarget <=0)
	{
		std::cerr << SUBPROGRAM ": invalid number of next target: " << opt::numOfNextTarget << ", must be greater than zero\n";
		die = true;
	}
	
	if(opt::maxLeaves <= 0)
	{
		std::cerr << SUBPROGRAM ":invalid number of max leaves:" << opt::maxLeaves << ", must be greater than zero\n";
		die = true;
	}
	
	if(opt::idmerLength <= 0)
	{
		std::cerr << SUBPROGRAM ":invalid kmer length to identify similar reads" << opt::idmerLength << ", must be greater than zero\n";
		die = true;
	}
	
	if(opt::minKmerLength <= 0)
	{
		std::cerr << SUBPROGRAM ":invalid min kmer length:" << opt::minKmerLength << ", must be greater than zero\n";
		die = true;
	}
	
	if(opt::genome != 5 && opt::genome != 10 && opt::genome != 100)
	{
		std::cerr << SUBPROGRAM ": invalid genome size: " << opt::genome << ", must be (5/10/100)[m]\n";
		die = true;
	}
	
	if(opt::mode != 0 && opt::mode != 1 && opt::mode != 2)
	{
		std::cerr << SUBPROGRAM ": invalid mode: " << opt::mode << ", must be (0/1/2)\n";
		die = true;
	}
	
	if(die)
	{
		std::cerr << "\n" << CORRECT_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}
	
	opt::readsFile = argv[optind++];
	//std::string out_prefix = stripFilename(opt::readsFile);
	//opt::correctFile = opt::directory + out_prefix + ".correct.fa";
	//opt::discardFile = opt::directory + out_prefix + ".discard.fa";
	opt::correctFile = opt::directory + "correct.fa";
	opt::discardFile = opt::directory + "discard.fa";
	
}
