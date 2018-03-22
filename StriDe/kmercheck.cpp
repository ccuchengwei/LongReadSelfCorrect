//
// kmercheck - get sequences kmer frequency
//
#include "SGACommon.h"
#include "Util.h"
#include "kmercheck.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "SequenceProcessFramework.h"
#include "KmerCheckProcess.h"
#include <iostream>
#include <map>
#include <memory>
//
// Getopt
//
#define SUBPROGRAM "kmercheck"

static const char *KMERFREQ_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"\n";

static const char *KMERFREQ_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Get sequences kmer frequency\n"
"  -p, --prefix=PREFIX       Use PREFIX for the names of the index files\n"
"  -o, --directory=PATH      Put results in the directory\n"
"  -t, --threads=NUM         Use NUM threads for the computation (default: 1)\n"
"  -r, --reference=SEQ       The reference file\n"
"  -a, --align=PATH          Get alignment files in the directory\n"
"  -v, --verbose             Display verbose output\n"
"      --help                Display this help and exit\n"
"      --version             Display version\n";

static const char* PROGRAM_IDENT =
PACKAGE_NAME "::" SUBPROGRAM;

namespace opt
{
    static unsigned int verbose;
    static std::string prefix;
	static std::string directory;
	static int numThreads = 1;
	static std::string reference;
	static std::string align;
	
	static int sizeLb = 15;
	static int sizeUb = 15;	
	static std::string readsFile;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
}

static const char* shortopts = "p:o:t:r:a:v";

enum { OPT_ALIGN, OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "prefix",      required_argument, nullptr, 'p' },
	{ "directory",   required_argument, nullptr, 'o' },
	{ "threads",     required_argument, nullptr, 't' },
	{ "reference",   required_argument, nullptr, 'r' },
	{ "align",       required_argument, nullptr, 'a' },
    { "verbose",     no_argument,       nullptr, 'v' },
    { "help",        no_argument,       nullptr, OPT_HELP },
    { "version",     no_argument,       nullptr, OPT_VERSION },
    { nullptr, 0, nullptr, 0 }
};

int kmercheckMain(int argc, char** argv)
{
	parseKMERCHECKOptions(argc, argv);
	
	BWT *pBWT, *pRBWT;
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
	}
	BWTIndexSet indexSet;
	indexSet.pBWT = pBWT;
	indexSet.pRBWT = pRBWT;
		
	//Set the kmer freq parameters
	//Load reference sequences
	std::unique_ptr<std::map<std::string, std::string> > pSeqMap(new std::map<std::string, std::string>);
	{
		SeqReader reader(opt::reference);
		WorkItemGenerator<SequenceWorkItem> generator(&reader);
		SequenceWorkItem workItem;
		while(generator.generate(workItem))
		{
			std::string id  = workItem.read.id;
			std::string seq = workItem.read.seq.toString();
			(*pSeqMap)[id] = seq;
		}
	}
	KmerCheckParameters kcParams;
	kcParams.indices     = indexSet;
	kcParams.directory   = opt::directory;
	kcParams.align       = opt::align;
	kcParams.pSeqMap     = pSeqMap.get();
	kcParams.size.first  = opt::sizeLb;
	kcParams.size.second = opt::sizeUb;
	
	std::cerr << "Using kmer size : " << kcParams.size.first << " - " << kcParams.size.second << "\n";

	Timer* pTimer = new Timer(PROGRAM_IDENT);
	
	KmerCheckPostProcess* pPostProcessor = new KmerCheckPostProcess(kcParams);
	
	if(opt::numThreads <= 1)
	{
		// Serial mode
		
		KmerCheckProcess* pProcessor = new KmerCheckProcess(kcParams);

		SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
		KmerCheckResult,
		KmerCheckProcess,
		KmerCheckPostProcess>(opt::readsFile, pProcessor, pPostProcessor);
		delete pProcessor;
		
	}
	else
	{
		// Parallel mode
		
		std::vector<KmerCheckProcess*> pProcessorVec;
		for(int i = 0; i < opt::numThreads; ++i)
			pProcessorVec.push_back(new KmerCheckProcess(kcParams));

		SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
		KmerCheckResult,
		KmerCheckProcess,
		KmerCheckPostProcess>(opt::readsFile, pProcessorVec, pPostProcessor);
		
		for(auto& iter : pProcessorVec)
			delete iter;
		
	}
	delete pPostProcessor;
	delete pBWT;
	delete pRBWT;
	delete pTimer;
	
    return 0;
}

// 
// Handle command line arguments
//
void parseKMERCHECKOptions(int argc, char** argv)
{
	optind = 1;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;) 
    {
        std::istringstream arg(optarg != nullptr ? optarg : "");
        switch (c) 
        {
			case 'p': arg >> opt::prefix; break;
			case 'o': arg >> opt::directory; break;
			case 't': arg >> opt::numThreads; break;
			case 'r': arg >> opt::reference; break;
			case 'a': arg >> opt::align; break;
			case 'v': opt::verbose++; break;
			case '?': die = true; break;
            case OPT_HELP:
                std::cout << KMERFREQ_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cout << KMERFREQ_VERSION_MESSAGE;
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
		if( system(("mkdir -p " + opt::directory + "split/").c_str()) != 0)
		{
			std::cerr << SUBPROGRAM << ": something wrong in directory: " << opt::directory << "\n";
			die = true;
		}
	}
	
	if(opt::numThreads <= 0)
	{
		std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::numThreads << "\n";
		die = true;
	}
	
	if(opt::reference.empty())
	{
		std::cerr << SUBPROGRAM << ": no reference\n";
		die = true;
	}
	
	if(opt::align.empty())
	{
		std::cerr << SUBPROGRAM << ": no alignments\n";
		die = true;
	}
	else
		opt::align += "/";
	
    if (die) 
    {
        std::cout << "\n" << KMERFREQ_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
	//Enter start and end kmersize
	while(true)
	{
		std::cout << "Please enter start and end kmer size\n";
		std::cin >> opt::sizeLb >> opt::sizeUb;
		if	(opt::sizeLb >= 10 && opt::sizeLb <= opt::sizeUb)
			break;
		std::cout << "Illegal values\n";
	}
	opt::readsFile = argv[optind++];
}
