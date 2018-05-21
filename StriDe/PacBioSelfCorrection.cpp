///-----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// PacBioSelfCorrectionProcess.cpp - Correction of PacBio reads using FM-index
//

#include <iostream>
#include <array>
#include <set>
#include "Util.h"
#include "PacBioSelfCorrection.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "SGACommon.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "SequenceProcessFramework.h"
#include "PacBioSelfCorrectionProcess.h"
#include "BWTIntervalCache.h"
#include "KmerThreshold.h"
#include "LongReadProbe.h"
#include "KmerCheckProcess.h"
#include "KmerFeature.h"
#include "BCode.h"

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
"      -t, --thread=NUM                 Use NUM threads for the computation (default: 1)\n"
"      -p, --prefix=PREFIX              Use PREFIX for the names of the index files\n"
"      -o, --output=DIR                 Output results in the directory\n"
"      -b, --barcode=FILE               Barcode of raw reads\n"
"\nPacBio correction parameters:\n"
"      -c, --PBcoverage=N               Coverage of PacBio reads (default: 90)\n"
"      -e, --error-rate=N               The error rate of PacBio reads.(default:0.15)\n"
"      -k, --kmer-size=N                The start kmer length (default: 19 (PacBioS).)\n"
"      -n, --next-target                The number of next FMWalk target seed(default: 1)\n"
"      -l, --max-leaves=N               Number of maximum leaves in the search tree. (default: 32)\n"
"      -i, --idmer-length=N             The length of the kmer to identify similar reads.(default: 9)\n"
"      -s, --min-kmer-size=N            The minimum length of the kmer to use. (default: 13.)\n"
"      -g, --genome=(5/10/100)[m]       Genome size of the species (default: 10m)\n"
"      -m, --mode=(0/1/2)               Mode in seed-searching (default: 1)\n"
"      -v, --verbose                    Display verbose output\n"
"      --help                           Display this help and exit\n"
"      --version                        Display version and exit\n"
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
	static int thread = 1;
	static std::string prefix;
	static std::string directory;
	static std::string barcode;
	static std::string readsFile;
	static std::string correctFile;
	static std::string discardFile;
	static BWTIndexSet indices;
	static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
	
    static size_t PBcoverage = 90; // PB seed searh depth
    static double ErrorRate = 0.15;
	
	static int startKmerLen = 19;
	static int nextTarget = 1;
	static int maxLeaves = 32;
    static int idmerLen = 9;
	static int minKmerLen = 13;
	
	static int genome = 10;
	static int mode = 1;
	static int verbose = 0;
	
	static bool Split = false;
    static bool DebugExtend = false;
    static bool DebugSeed = false;
	static bool OnlySeed = false;
	static bool NoDp = false;
	static bool Manual = false;
	
	//variables for auto set
	static bool Adjust = false;
	static std::map<int, int> order = { {5, 0}, {10, 1}, {100, 2} };
	static int size[3] = { 17, 19, 21 };
	static std::array<int, 3> offset = { 0, 0, 0 };
	static std::set<int> pool = { 5, 9, 19 }; // (5 & 9) -> fragment; 19 -> scan-window
}

static const char* shortopts = "t:p:o:b:c:e:k:u:r:n:l:i:s:g:m:v";

enum { OPT_HELP = 1, OPT_VERSION, OPT_SPLIT, OPT_FIRST, OPT_DEBUGEXTEND, OPT_DEBUGSEED, OPT_ONLYSEED, OPT_NODP };

static const struct option longopts[] = {
	{ "thread",             required_argument, nullptr, 't' },
	{ "prefix",             required_argument, nullptr, 'p' },
	{ "output",             required_argument, nullptr, 'o' },
	{ "barcode",            required_argument, nullptr, 'b' },
    { "PBcoverage",         required_argument, nullptr, 'c' },
    { "error-rate",         required_argument, nullptr, 'e' },
	{ "kmer-size",          required_argument, nullptr, 'k' },
	{ "unique-offset",      required_argument, nullptr, 'u' },
	{ "repeat-offset",      required_argument, nullptr, 'r' },
	{ "next-target",        required_argument, nullptr, 'n' },
	{ "max-leaves",         required_argument, nullptr, 'l' },
    { "idmer-length",       required_argument, nullptr, 'i' },
	{ "min-kmer-size",      required_argument, nullptr, 's' },
    { "genome",             required_argument, nullptr, 'g' },
    { "mode",               required_argument, nullptr, 'm' },
	{ "verbose",            no_argument,       nullptr, 'v' },
	{ "help",               no_argument,       nullptr, OPT_HELP },
	{ "version",            no_argument,       nullptr, OPT_VERSION },
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
	
	// Load indices
	std::unique_ptr<BWT> pBWT, pRBWT;
	std::unique_ptr<SampledSuffixArray> pSSA;
	#pragma omp parallel sections
	{
		#pragma omp section
		{	//Initialization of large BWT takes some time, pass the disk to next job
			std::cerr << "Loading BWT: " << opt::prefix + BWT_EXT << "\n";
			pBWT = std::unique_ptr<BWT>(new BWT(opt::prefix + BWT_EXT, opt::sampleRate));
		}
		#pragma omp section
		{
			std::cerr << "Loading RBWT: " << opt::prefix + RBWT_EXT << "\n";
			pRBWT = std::unique_ptr<BWT>(new BWT(opt::prefix + RBWT_EXT, opt::sampleRate));
		}
		#pragma omp section
		{
			std::cerr << "Loading Sampled Suffix Array: " << opt::prefix + SAI_EXT << "\n";
			pSSA = std::unique_ptr<SampledSuffixArray>(new SampledSuffixArray(opt::prefix + SAI_EXT, SSA_FT_SAI));
		}
	}
	opt::indices.pBWT  = pBWT.get();
	opt::indices.pRBWT = pRBWT.get();
	opt::indices.pSSA  = pSSA.get();
	
	ecParams.indices   = opt::indices;
	ecParams.directory = opt::directory;

    ecParams.PBcoverage = opt::PBcoverage;
    ecParams.ErrorRate  = opt::ErrorRate;
	ecParams.nextTarget = opt::nextTarget;
	ecParams.maxLeaves  = opt::maxLeaves;
    ecParams.idmerLen   = opt::idmerLen;
	ecParams.minKmerLen = opt::minKmerLen;
	
	ecParams.Split       = opt::Split;
    ecParams.DebugExtend = opt::DebugExtend;
    ecParams.DebugSeed   = opt::DebugSeed;
	
	if(opt::OnlySeed) BCode::load(opt::barcode);
	ecParams.OnlySeed    = opt::OnlySeed;
	ecParams.NoDp        = opt::NoDp;
	
	if(!opt::Adjust)
	{
		opt::startKmerLen  = opt::size[opt::order[opt::genome]];
		opt::offset[1] = 2 * std::min(std::max((ecParams.PBcoverage/30 - 1), 0), (opt::order[opt::genome] + 1));
		opt::offset[2] = -2 * (opt::order[opt::genome] + 1);
	}
	ecParams.startKmerLen = opt::startKmerLen;
	
	//Insert kmer sizes in pool for future usage. Noted by KuanWeiLee 18/3/12
	for(auto& o : opt::offset)
		opt::pool.insert(opt::startKmerLen + o);
	ecParams.pool = opt::pool;
	
	FMextendParameters FM_params(
			opt::indices,
			opt::idmerLen,
			opt::maxLeaves,
			opt::minKmerLen,
			opt::PBcoverage,
			opt::ErrorRate);
	ecParams.FM_params = FM_params;
	
	
	LongReadProbe::m_params = 
	ProbeParameters(
			opt::indices,
			opt::directory,
			opt::startKmerLen,
			opt::PBcoverage,
			opt::mode,
			opt::offset,
			opt::pool,
			opt::DebugSeed,
			opt::Manual);
	
	//Initialize KmerThreshold
	KmerThreshold::Instance().initialize(-1, 50, opt::PBcoverage, opt::directory);
	
	std::cerr
	<< "\nCorrecting PacBio reads for " << opt::readsFile << " using--\n"
	<< "number of threads:\t" << opt::thread << "\n"
	<< "PB reads coverage:\t" << opt::PBcoverage << "\n"
	<< "num of next Targets:\t" << opt::nextTarget << "\n"
	<< "large kmer size:\t" << opt::startKmerLen << "\n" 
	<< "small kmer size:\t" << opt::minKmerLen << "\n"
	<< "max leaves:\t" << opt::maxLeaves  << "\n"
	<< "max depth:\t1.2~0.8* (length between two seeds +- 20)" << "\n";

	
	// Start a timer
	Timer* pTimer = new Timer(PROGRAM_IDENT);
	
	//Start processing sequences
	SequenceProcessFramework::processSequences<SequenceWorkItem,
	PacBioSelfCorrectionResult,
	PacBioSelfCorrectionProcess,
	PacBioSelfCorrectionPostProcess,
	PacBioSelfCorrectionParameters>(opt::thread, opt::readsFile, ecParams);
	
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
			case 't': arg >> opt::thread;     break;
			case 'p': arg >> opt::prefix;     break;
			case 'o': arg >> opt::directory;  break;
			case 'b': arg >> opt::barcode;    break;
			case 'c': arg >> opt::PBcoverage; break;
			case 'e': arg >> opt::ErrorRate;  break;
			case 'k':
				arg >> opt::startKmerLen;
				opt::Adjust = true;
				break;
			case 'u':
			//	arg >> opt::unique_offset;
				arg >> opt::offset[1];
				opt::Adjust = true;
				break;
			case 'r':
			//	arg >> opt::repeat_offset;
				arg >> opt::offset[2];
				opt::Adjust = true;
				break;
			case 'n': arg >> opt::nextTarget; break;
			case 'l': arg >> opt::maxLeaves;  break;
			case 'i': arg >> opt::idmerLen;   break;
			case 's': arg >> opt::minKmerLen; break;
			case 'g': arg >> opt::genome;     break;
			case 'm':
				arg >> opt::mode;
				opt::Manual = true;
				break;
			case 'v': opt::verbose++; break;
			case OPT_HELP:
				std::cerr << CORRECT_USAGE_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_VERSION:
				std::cerr << CORRECT_VERSION_MESSAGE;
				exit(EXIT_SUCCESS);
			case OPT_SPLIT:       opt::Split       = true; break;
			case OPT_DEBUGEXTEND: opt::DebugExtend = true; break;
			case OPT_DEBUGSEED:   opt::DebugSeed   = true; break;
			case OPT_NODP:        opt::NoDp        = true; break;
			case OPT_ONLYSEED:
				opt::DebugSeed = true;
				opt::OnlySeed = true;
				break;
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

	if(opt::thread <= 0)
	{
		std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::thread << "\n";
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
		std::vector<std::string> subdir;
		subdir.push_back(std::string(""));
		if(opt::DebugSeed)
		{
			subdir.pop_back();
			subdir.push_back(std::string("extend/"));
			subdir.push_back(std::string("seed/error/"));
		}
		for(auto& iter : subdir)
			if(system(("mkdir -p " + opt::directory + iter).c_str()) != 0)
			{
				std::cerr << SUBPROGRAM << ": something wrong making directory: " << opt::directory << "\n";
				die = true;
			}
	}

	if(opt::PBcoverage <= 0)
	{
		std::cerr << SUBPROGRAM ": invalid number of coverage: " << opt::PBcoverage << ", must be greater than zero\n";
		die = true;
	}
	
	if(opt::ErrorRate < 0 || opt::ErrorRate > 1)
	{
		std::cerr <<SUBPROGRAM ":invalid error rate: " << opt::ErrorRate << ", must be 0 ~ 1\n";
		die = true;
	}
	
	if(opt::startKmerLen <= 0)
	{
		std::cerr << SUBPROGRAM ": invalid start kmer length: " << opt::startKmerLen << ", must be greater than zero\n";
		die = true;
	}
	
	if(opt::nextTarget <= 0)
	{
		std::cerr << SUBPROGRAM ": invalid number of next target: " << opt::nextTarget << ", must be greater than zero\n";
		die = true;
	}
	
	if(opt::maxLeaves <= 0)
	{
		std::cerr << SUBPROGRAM ":invalid number of max leaves:" << opt::maxLeaves << ", must be greater than zero\n";
		die = true;
	}
	
	if(opt::idmerLen <= 0)
	{
		std::cerr << SUBPROGRAM ":invalid kmer length to identify similar reads" << opt::idmerLen << ", must be greater than zero\n";
		die = true;
	}
	
	if(opt::minKmerLen <= 0)
	{
		std::cerr << SUBPROGRAM ":invalid min kmer length:" << opt::minKmerLen << ", must be greater than zero\n";
		die = true;
	}
	
	if(opt::genome != 5 && opt::genome != 10 && opt::genome != 100)
	{
		std::cerr << SUBPROGRAM ": invalid genome size: " << opt::genome << ", must be (5/10/100)[m]\n";
		die = true;
	}
	
	if(opt::mode < 0 || opt::mode > 2)
	{
		std::cerr << SUBPROGRAM ": invalid mode: " << opt::mode << ", must be (0/1/2)\n";
		die = true;
	}
	
	if(opt::OnlySeed && opt::barcode.empty())
	{
		std::cerr << SUBPROGRAM ": no barcode\n";
		die = true;
	}
	
	if(die)
	{
		std::cerr << "\n" << CORRECT_USAGE_MESSAGE;
		exit(EXIT_FAILURE);
	}
	
	opt::readsFile = argv[optind++];
//	std::string out_prefix = stripFilename(opt::readsFile);
	
}

