//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// KmerDistribution - Histogram of kmer frequencies
//
#include "KmerDistribution.h"
#include <limits>
#include <cmath>
#include <fstream>
void KmerDistribution::operator+=(const KmerDistribution& temp)
{
	for(iteratorKmerFreqsMap iter = temp.m_data.begin(); iter != temp.m_data.end(); iter++)
		m_data[iter->first] += iter->second;
}
double KmerDistribution::getCumulativeProportionLEQ(int n) const
{    
	size_t cumulativeSum = 0;
	double cumulativeProportion = 0;
	for(iteratorKmerFreqsMap iter = m_data.begin(); iter != m_data.end(); iter++)
	{
		if(iter->first > n)
			break;
		cumulativeSum +=  iter->second;
		cumulativeProportion = (double)cumulativeSum/m_total;
	}
	return cumulativeProportion;
}

size_t KmerDistribution::getCutoffForProportion(double p) const
{
	if(p > 1 || p < 0)
	{
		std::cout<<"Portion should between 0 <-> 1.\n";
		exit(EXIT_FAILURE);
	}
	size_t kmerFreqs = 0;
	size_t cumulativeSum = 0;
	double cumulativeProportion = 0;
	for(iteratorKmerFreqsMap iter = m_data.begin(); iter != m_data.end(); iter++)
	{
		kmerFreqs = iter->first;
		cumulativeSum += iter->second;
		cumulativeProportion = (double)cumulativeSum/m_total;
		if(cumulativeProportion > p)
			break;
	}
	return kmerFreqs;
}
//compute median and std
void KmerDistribution::computeKDAttributes(float censor)
{
	std::vector<size_t>rawdata;
	iteratorKmerFreqsMap iter = m_data.begin();
	while(iter != m_data.end() && iter->first <= censor)
		iter++;
	for(size_t most = 0; iter !=  m_data.end(); iter++)
	{
		rawdata.resize((rawdata.size() + iter->second), iter->first);
		if(iter->second > most)
		{
			m_mode = iter->first;
			most = iter->second;
		}
	}
    
	//compute quartiles
	if(rawdata.empty())return;
	m_median = rawdata[rawdata.size()/2];
	m_first_quartile = rawdata[rawdata.size()/4];
	m_third_quartile = rawdata[rawdata.size()*3/4];
	
	//compute standard deviation
	if(rawdata.size() <= 1)return;
	double difference_square = 0;
	for(std::vector<size_t>::iterator iterFreq = rawdata.begin(); iterFreq != rawdata.end(); iterFreq++)
		difference_square += pow(((*iterFreq) - m_median), 2);
	double variance = difference_square / (rawdata.size() - 1);
	m_std = sqrt(variance);
	
	// double freq95 = getCutoffForProportion(0.95);
	// m_repeatKmerCutoff = m_std > m_median*2? m_median*1.5: m_median*1.3;
	m_repeatKmerCutoff = m_median*1.3;
	// m_repeatKmerCutoff = getCutoffForProportion(0.8);
	// m_repeatKmerCutoff = (double) m_median*(0.39+0.53* (freq95/(double)m_median));
	
}
void KmerDistribution::write(std::ostream& outfile,TYPE mode) const
{
	switch(mode)
	{
		case DATA:
			for(iteratorKmerFreqsMap iter = m_data.begin(); iter !=  m_data.end(); iter++)
				outfile << iter->first << "\t" << iter->second << "\n";
			break;
		case ATTRIBUTE:
			outfile << m_readid << "\t " << m_first_quartile << "\t" << m_median << "\t" << m_third_quartile << "\t" << m_mode << "\t" << m_std <<"\n";
			break;
		case SPLIT:
			for(iteratorKmerFreqsMap iter = m_data.begin(); iter != m_data.end(); iter++)
				for(size_t i = 1; i <= iter->second; i++)
					outfile << iter->first << "\n";
			break;
		default:
			std::cout << "KmerDistribution-write-mode : DATA , ATTRIBUTE\n";
			exit(EXIT_FAILURE);
	}
}

//Below are legacy codes. Noted by KuanWeiLee 20171027
/***********************************************************************************/
int KmerDistribution::findFirstLocalMinimum() const
{
    std::vector<int> countVector = toCountVector(1000);
    if(countVector.empty())
        return -1;

    std::cout << "CV: " << countVector.size() << "\n";
    int prevCount = countVector[1];
    double threshold = 0.75;
    for(size_t i = 2; i < countVector.size(); ++i)
    {
        int currCount = countVector[i];
        double ratio = (double)currCount / prevCount;
        std::cout << i << " " << currCount << " " << ratio << "\n";
        if(ratio > threshold)
            return i;
        prevCount = currCount;
    }
    return -1;
}


// Find the boundary of the kmers that are likely erroneous
// We do this by finding the value between 1 and the trusted mode
// that contributes the fewest
int KmerDistribution::findErrorBoundary() const
{
    int mode = getCensoredMode(5);
    if(mode == -1)
        return -1;

    std::cerr << "Trusted kmer mode: " << mode  << "\n";
    std::vector<int> countVector = toCountVector(1000);
    if(countVector.empty())
        return -1;

    int runningSum = 0;
    double minContrib = std::numeric_limits<double>::max();
    int idx = -1;
    for(int i = 1; i < mode; ++i)
    {
        runningSum += countVector[i];
        double v = (double)countVector[i] / runningSum;
        if(v < minContrib)
        {
            minContrib = v;
            idx = i;
        }
    }
    return idx;
}

//Similar to findFirstLocalMinimum() with a few different parameters 
int KmerDistribution::findErrorBoundaryByRatio(double ratio) const
{
    int mode = getCensoredMode(5);
    if(mode == -1)
        return -1;

    std::cerr << "Trusted kmer mode: " << mode  << "\n";
    std::vector<int> countVector = toCountVector(1000);
    if(countVector.empty())
        return -1;

    for(int i = 1; i < mode - 1; ++i)
    {
        int currCount = countVector[i];
        int nextCount  = countVector[i+1];
        double cr = (double)currCount / nextCount;
        if(cr < ratio)
            return i;
    }
    return -1;
}
size_t KmerDistribution::getCensoredMode(size_t n) const
{
	size_t most = 0,mode = 0;
    iteratorKmerFreqsMap iter = m_data.begin();
	while(iter != m_data.end() && iter->first < n)
		iter++;
	for(; iter != m_data.end(); iter++)
		if(iter->second > most)
		{
			mode = iter->first;
			most = iter->second;
		}
	return mode;
}
// 
std::vector<int> KmerDistribution::toCountVector(int max) const
{
    std::vector<int> out;
    if(m_data.empty())
        return out;

    int min = 0;

    for(int i = min; i <= max; ++i)
    {
        std::map<size_t,size_t>::const_iterator iter = m_data.find(i);
        int v = (iter != m_data.end()) ? iter->second : 0;
        out.push_back(v);
    }
    return out;
}

// for compatibility with old code
void KmerDistribution::print(int max) const
{
    print(stdout, max);
}

void KmerDistribution::print(FILE* fp, int max) const
{
    fprintf(fp, "Kmer coverage histogram\n");
    fprintf(fp, "cov\tcount\n");

    int maxCount = 0;
    std::map<size_t,size_t>::const_iterator iter = m_data.begin();
    for(; iter != m_data.end(); ++iter)
    {
        if(iter->first <= max)
            fprintf(fp, "%d\t%d\n", iter->first, iter->second);
        else
            maxCount += iter->second;
    }
    fprintf(fp, ">%d\t%d\n", max, maxCount);

}

