//
#include "SGACommon.h"
#include "Util.h"
#include "checkkmer.h"
#include "SuffixArray.h"
#include "BWT.h"
#include "Timer.h"
#include "BWTAlgorithms.h"
#include "SequenceProcessFramework.h"
#include "CheckKmerProcess.h"
#include <iostream>
#include <memory>
#include <map>
#include <list>
//
// Getopt
//
#define SUBPROGRAM "checkkmer"

static const char *KMERFREQ_VERSION_MESSAGE =
SUBPROGRAM " Version " PACKAGE_VERSION "\n"
"\n";

static const char *KMERFREQ_USAGE_MESSAGE =
"Usage: " PACKAGE_NAME " " SUBPROGRAM " [OPTION] ... READSFILE\n"
"Get sequences kmer frequency\n"
"  -p, --prefix=PREFIX       Use PREFIX for the names of the index files\n"
"  -o, --directory=PATH      Put results in the directory\n"
"  -t, --threads=NUM         Use NUM threads for the computation (default: 1)\n"
"  -b, --barcode=FILE        Use the barcode to check kmer \n"
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
	static int thread = 1;
	static std::string barcode;
	static int sizeLb = 15;
	static int sizeUb = 15;	
	static std::string readsFile;
    static int sampleRate = BWT::DEFAULT_SAMPLE_RATE_SMALL;
}

static const char* shortopts = "p:o:t:b:v";

enum { OPT_ALIGN, OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
	{ "prefix",      required_argument, nullptr, 'p' },
	{ "directory",   required_argument, nullptr, 'o' },
	{ "threads",     required_argument, nullptr, 't' },
	{ "barcode",     required_argument, nullptr, 'b' },
    { "verbose",     no_argument,       nullptr, 'v' },
    { "help",        no_argument,       nullptr, OPT_HELP },
    { "version",     no_argument,       nullptr, OPT_VERSION },
    { nullptr, 0, nullptr, 0 }
};

int checkkmerMain(int argc, char** argv)
{
	parseCHECKKMEROptions(argc, argv);
	
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
	
	std::unique_ptr<std::map<std::string, std::list<CodeBlock> > > pAlignLog(new std::map<std::string, std::list<CodeBlock> >);
	std::cerr << "Loading BARCODE: " << opt::barcode << '\n';
	std::istream* pCodeReader = createReader(opt::barcode);
	while(true)
	{
		if(pCodeReader->eof()) break;
		std::string qname, tname, code, rvc, sup;
		int qstart, qend, tstart, tend;
		*pCodeReader
		>> qname >> qstart >> qend
		>> tname >> tstart >> tend
		>> code  >> rvc    >> sup;
		(*pAlignLog)[qname].push_back(CodeBlock(qstart, qend, code, (rvc == "True" ? true : false)));
	}
	delete pCodeReader;
	
	BWTIndexSet indexSet;
	indexSet.pBWT  = pBWT.get();
	indexSet.pRBWT = pRBWT.get();
	
	
	CheckKmerParameters ckParams;
	ckParams.indices     = indexSet;
	ckParams.directory   = opt::directory;
	ckParams.size.first  = opt::sizeLb;
	ckParams.size.second = opt::sizeUb;
	ckParams.pAlignLog   = pAlignLog.get();
	ckParams.mode        = true;
	
	std::cerr << "Using kmer size : " << ckParams.size.first << " - " << ckParams.size.second << "\n";
	
	Timer* pTimer = new Timer(PROGRAM_IDENT);
	
	SequenceProcessFramework::processSequences<SequenceWorkItem,
	CheckKmerResult,
	CheckKmerProcess,
	CheckKmerPostProcess,
	CheckKmerParameters>(opt::thread, opt::readsFile, ckParams);
	
	delete pTimer;
    return 0;
}

// 
// Handle command line arguments
//
void parseCHECKKMEROptions(int argc, char** argv)
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
			case 't': arg >> opt::thread; break;
			case 'b': arg >> opt::barcode; break;
			case 'v': opt::verbose++; break;
			case '?': die = true; break;
            case OPT_HELP:
                std::cerr << KMERFREQ_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cerr << KMERFREQ_VERSION_MESSAGE;
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
		if(system(("mkdir -p " + opt::directory + "split/").c_str()) != 0)
		{
			std::cerr << SUBPROGRAM << ": something wrong in directory: " << opt::directory << "\n";
			die = true;
		}
	}
	
	if(opt::thread <= 0)
	{
		std::cerr << SUBPROGRAM ": invalid number of threads: " << opt::thread << "\n";
		die = true;
	}
	
	if(opt::barcode.empty())
	{
		std::cerr << SUBPROGRAM << ": no barcode\n";
		die = true;
	}
	
    if (die) 
    {
        std::cerr << "\n" << KMERFREQ_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
	
	while(true)
	{
		std::cerr << "Please enter start & end kmer size\n";
		std::cin >> opt::sizeLb >> opt::sizeUb;
		if(opt::sizeLb >= 7 && opt::sizeUb >= opt::sizeLb) break;
		std::cerr << "Illegal values\n";
	}
	
	opt::readsFile = argv[optind++];
}
