#include "Kmer.h"

//Kmer Map for different size of kmers. Noted by KuanWeiLee 20180208
thread_local std::map<int, std::unique_ptr<Kmer[]> > Kmer::kmerMap;