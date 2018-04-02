//
// kmerfreq - get sequences kmer frequency
//
#include <iostream>
#include <memory>
#include "SGACommon.h"
#include "Util.h"
#include "kmerfreq.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "BWTAlgorithms.h"
#include "KmerFeature.h"
#include "KmerThreshold.h"
//
// Getopt
//
#define SUBPROGRAM "kmerfreq"

static const char *KMERFREQ_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"\n";

static const char *KMERFREQ_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION]\n"
"Get sequences kmer frequency\n"
"  -p, --prefix=PREFIX       Use PREFIX for the names of the index files\n"
"  -c, --PBcoverage=N        Coverage of PacBio reads (default: 90)\n"
"  -v, --verbose             Display verbose output\n"
"      --help                Display this help and exit\n"
"      --version             Display version\n";

namespace opt
{
    static unsigned int verbose;
	static std::string prefix;
	static int PBcoverage = 90;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
}

static const char* shortopts = "p:c:v";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "prefix",      required_argument, nullptr, 'p' },
	{ "PBcoverage",  required_argument, nullptr, 'c' },
    { "verbose",     no_argument,       nullptr, 'v' },
    { "help",        no_argument,       nullptr, OPT_HELP },
    { "version",     no_argument,       nullptr, OPT_VERSION },
    { nullptr, 0, nullptr, 0 }
};

int kmerfreqMain(int argc, char** argv)
{
	parseKMERFREQOptions(argc, argv);
	
	std::unique_ptr<BWT> pBWT, pRBWT;
	// Load indices
	#pragma omp parallel
	{
		#pragma omp single nowait
		{	//Initialization of large BWT takes some time, pass the disk to next job
			std::cerr << "Loading BWT: " << opt::prefix + BWT_EXT << "\n";
			pBWT = std::unique_ptr<BWT>(new BWT(opt::prefix + BWT_EXT, opt::sampleRate));
		}
		#pragma omp single nowait
		{
			std::cerr << "Loading RBWT: " << opt::prefix + RBWT_EXT << "\n";
			pRBWT = std::unique_ptr<BWT>(new BWT(opt::prefix + RBWT_EXT, opt::sampleRate));
		}
	}
	BWTIndexSet indexSet;
	indexSet.pBWT = pBWT.get();
	indexSet.pRBWT = pRBWT.get();
	KmerThreshold::Instance().set(-1, 100, opt::PBcoverage, "");
	
	std::string query;
	int staticSize;
	int mode;
	std::cerr << "Please enter query sequence, kmer size and mode:\n";
	while (std::cin >> query >> staticSize >> mode) 
	{
		KmerFeature* prev = nullptr;
		int queryLen = query.length();
		int dynamicSize = staticSize;
		for(int pos = 0; pos <= (queryLen - staticSize); pos++)
		{
			KmerFeature staticKmer(indexSet, query, pos, staticSize);
			KmerFeature* dynamicKmer = new KmerFeature(indexSet, query, 0, dynamicSize++, prev);
			std::cout << pos << '\t'
			<< staticKmer.getWord() << '\t' << staticKmer.getFreq()
			<< " <-> " << KmerThreshold::Instance().get(mode, staticSize) << '\t'
			<< dynamicKmer->getWord() << '\t' << dynamicKmer->getFreq()
			<< " <-> " << KmerThreshold::Instance().get(mode, dynamicKmer->getSize()) << '\n';
			delete prev;
			prev = dynamicKmer;
		}
		std::cout << "-\n";
		delete prev;
	}
	std::cerr << "Exit successfully!\n";
    return 0;
}

// 
// Handle command line arguments
//
void parseKMERFREQOptions(int argc, char** argv)
{
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;) 
    {
        std::istringstream arg(optarg != nullptr ? optarg : "");
        switch (c) 
        {
			case 'p': arg >> opt::prefix; break;
			case 'c': arg >> opt::PBcoverage; break;
            case '?': die = true; break;
            case 'v': opt::verbose++; break;
            case OPT_HELP:
                std::cout << KMERFREQ_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << KMERFREQ_VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
	
	if(opt::prefix.empty())
	{
		std::cerr << SUBPROGRAM << ": no prefix\n";
		die = true;
	}
	
	if(opt::PBcoverage <= 0)
	{
		std::cerr << SUBPROGRAM ": invalid number of coverage: " << opt::PBcoverage << ", must be greater than zero\n";
		die = true;
	}
	
    if (die) 
    {
        std::cout << "\n" << KMERFREQ_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }

}
