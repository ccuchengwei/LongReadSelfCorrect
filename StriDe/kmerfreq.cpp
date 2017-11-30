//
// kmerfreq - get sequences kmer frequency
//
#include <iostream>
#include <map>
#include "SGACommon.h"
#include "Util.h"
#include "kmerfreq.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "SequenceProcessFramework.h"
#include "KmerFreqProcess.h"
//
// Getopt
//
#define SUBPROGRAM "kmerfreq"

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
	
	static int kmerSizeLb = 15;
	static int kmerSizeUb = 15;	
	static std::string readsFile;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
}

static const char* shortopts = "p:o:t:r:a:v";

enum { OPT_ALIGN, OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "prefix",      required_argument, NULL, 'p' },
	{ "directory",   required_argument, NULL, 'o' },
	{ "threads",     required_argument, NULL, 't' },
	{ "reference",   required_argument, NULL, 'r' },
	{ "align",       required_argument, NULL, 'a' },
    { "verbose",     no_argument,       NULL, 'v' },
    { "help",        no_argument,       NULL, OPT_HELP },
    { "version",     no_argument,       NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

int kmerfreqMain(int argc, char** argv)
{
	parseKMERFREQOptions(argc, argv);
	
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
	std::map<std::string,std::string>refMap;
	{
		SeqReader reader(opt::reference);
		WorkItemGenerator<SequenceWorkItem> generator(&reader);
		SequenceWorkItem workItem;
		std::string id,seq;
		while(generator.generate(workItem))
		{
			id = workItem.read.id;
			seq = workItem.read.seq.toString();
			refMap[id] = seq;
		}
	}
	KmerFreqParameters kfParams(refMap);
	kfParams.indices = indexSet;
	kfParams.directory = opt::directory;
	kfParams.align = opt::align;
	kfParams.kmerSize.first = opt::kmerSizeLb;
	kfParams.kmerSize.second = opt::kmerSizeUb;

	KmerFreqPostProcess* pPostProcessor = new KmerFreqPostProcess(kfParams);

	
	Timer* pTimer = new Timer(PROGRAM_IDENT);
	
	if(opt::numThreads <= 1)
	{
		// Serial mode
		
		KmerFreqProcess* pProcessor = new KmerFreqProcess(kfParams);

		SequenceProcessFramework::processSequencesSerial<SequenceWorkItem,
		KmerFreqResult,
		KmerFreqProcess,
		KmerFreqPostProcess>(opt::readsFile, pProcessor, pPostProcessor);
		delete pProcessor;
		
	}
	else
	{
		// Parallel mode
		
		std::vector<KmerFreqProcess*> pProcessorVector;
		for(int i = 0; i < opt::numThreads; ++i)
		{
			KmerFreqProcess* pProcessor = new KmerFreqProcess(kfParams);
			pProcessorVector.push_back(pProcessor);
		}

		SequenceProcessFramework::processSequencesParallel<SequenceWorkItem,
		KmerFreqResult,
		KmerFreqProcess,
		KmerFreqPostProcess>(opt::readsFile, pProcessorVector, pPostProcessor);
		
		while(!pProcessorVector.empty())
		{
			delete pProcessorVector.back();
			pProcessorVector.pop_back();
		}
		
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
void parseKMERFREQOptions(int argc, char** argv)
{
	optind = 1;
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) 
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
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
		std::cin >> opt::kmerSizeLb >> opt::kmerSizeUb;
		if	(opt::kmerSizeLb >= 10 && opt::kmerSizeLb <= opt::kmerSizeUb)
			break;
		std::cout << "Illegal values\n";
	}
	opt::readsFile = argv[optind++];
}
