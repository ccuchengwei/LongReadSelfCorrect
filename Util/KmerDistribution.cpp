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

void KmerDistribution::operator+=(const KmerDistribution& other)
{
	for(const auto& iter : other.data)
		this->data[iter.first] += iter.second;
}

int KmerDistribution::getNumberWithCount (int n) const
{
	std::map<int, int>::const_iterator iter = data.find(n);
	return iter != data.end() ? iter->second : 0;
}

int KmerDistribution::getQuartile(int n) const
{
	switch(n)
	{
		case 1: return first_quartile;
		case 2: return median;
		case 3: return third_quartile;
		default: 
			std::cout << "Quartile options : 1 , 2 ,3\n";
			exit(EXIT_FAILURE);
	}
}
		
double KmerDistribution::getCumulativeProportionLEQ(int n) const
{    
	int cumulativeSum = 0;
	double cumulativeProportion = 0;
	for(const auto& iter : data)
	{
		if(iter.first > n) break;
		cumulativeSum += iter.second;
		cumulativeProportion = (double)cumulativeSum/total;
	}
	return cumulativeProportion;
}

int KmerDistribution::getCutoffForProportion(double p) const
{
	if(p > 1 || p < 0)
	{
		std::cerr<<"Portion should between 0 <-> 1.\n";
		exit(EXIT_FAILURE);
	}
	int kmerFreq = 0;
	int cumulativeSum = 0;
	double cumulativeProportion = 0;
	for(const auto& iter : data)
	{
		kmerFreq = iter.first;
		cumulativeSum += iter.second;
		cumulativeProportion = (double)cumulativeSum/total;
		if(cumulativeProportion > p) break;
	}
	return kmerFreq;
}

//compute median and sdv
void KmerDistribution::computeKDAttributes(float censor)
{
	std::vector<int>rawdata;
	std::map<int, int>::const_iterator iter = data.begin();
	
	for(; iter != data.end() && iter->first <= censor; iter++);
	
	for(int most = 0; iter !=  data.end(); iter++)
	{
		rawdata.resize((rawdata.size() + iter->second), iter->first);
		if(iter->second > most)
		{
			mode = iter->first;
			most = iter->second;
		}
	}
	
	//compute quartiles
	if(rawdata.empty()) return;
	median = rawdata[rawdata.size()/2];
	first_quartile = rawdata[rawdata.size()/4];
	third_quartile = rawdata[rawdata.size()*3/4];
	
	//compute standard deviation
	if(rawdata.size() <= 1) return;
	double difference_square = 0;
	for(const auto& iter : rawdata)
		difference_square += pow((iter - median), 2);
	double variance = difference_square/(rawdata.size() - 1);
	sdv = sqrt(variance);
	
	// double freq95 = getCutoffForProportion(0.95);
	// repeatKmerCutoff = sdv > median*2? median*1.5: median*1.3;
	repeatKmerCutoff = median*1.3;
	// repeatKmerCutoff = getCutoffForProportion(0.8);
	// repeatKmerCutoff = (double) median*(0.39+0.53* (freq95/(double)median));
	
}

void KmerDistribution::write(std::ostream& out, int mode) const
{
	switch(mode)
	{
		case 0:
			for(const auto& iter : data)
				out << iter.first <<'\t' << iter.second << '\n';
			break;
		case 1:
			out << readid << '\t' << first_quartile << '\t' << median << '\t' << third_quartile << '\t' << mode << '\t' << sdv << '\n';
			break;
		default:
			std::cerr << "Mode should be 0/1.\n";
			exit(EXIT_FAILURE);
	}
}

//Legacy Part
/***********/
int KmerDistribution::findFirstLocalMinimum() const
{
	std::vector<int> countVector = toCountVector(1000);
	if(countVector.empty()) return -1;

	std::cout << "CV: " << countVector.size() << '\n';
	int prevCount = countVector[1];
	double threshold = 0.75;
	for(int i = 2; i < (int)countVector.size(); ++i)
	{
		int currCount = countVector[i];
		double ratio = (double)currCount / prevCount;
		std::cout << i << " " << currCount << " " << ratio << '\n';
		if(ratio > threshold) return i;
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
	if(mode == -1) return -1;

	std::cerr << "Trusted kmer mode: " << mode  << '\n';
	std::vector<int> countVector = toCountVector(1000);
	if(countVector.empty()) return -1;

	int runningSum = 0;
	double minContrib = std::numeric_limits<double>::max();
	int idx = -1;
	for(int i = 1; i < mode; ++i)
	{
		runningSum += countVector[i];
		double v = (double)countVector[i]/runningSum;
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
	if(mode == -1) return -1;

	std::cerr << "Trusted kmer mode: " << mode  << '\n';
	std::vector<int> countVector = toCountVector(1000);
	if(countVector.empty())
		return -1;

	for(int i = 1; i < mode - 1; ++i)
	{
		int currCount = countVector[i];
		int nextCount  = countVector[i+1];
		double cr = (double)currCount / nextCount;
		if(cr < ratio) return i;
	}
	return -1;
}

int KmerDistribution::getCensoredMode(int n) const
{
	int most = 0, mode = 0;
	std::map<int, int>::const_iterator iter = data.begin();
	
	for(; iter != data.end() && iter->first < n; iter++);
	
	for(; iter != data.end(); iter++)
		if(iter->second > most)
		{
			mode = iter->first;
			most = iter->second;
		}
	return mode;
}

std::vector<int> KmerDistribution::toCountVector(int max) const
{
	std::vector<int> out;
	if(data.empty())
		return out;

	int min = 0;

	for(int i = min; i <= max; ++i)
	{
		std::map<int, int>::const_iterator iter = data.find(i);
		int v = iter != data.end() ? iter->second : 0;
		out.push_back(v);
	}
	return out;
}

void KmerDistribution::print(int max) const
{
	print(stdout, max);
}

void KmerDistribution::print(FILE *fp, int max) const
{
	fprintf(fp, "Kmer coverage histogram\n");
	fprintf(fp, "cov\tcount\n");

	int maxCount = 0;
	std::map<int, int>::const_iterator iter = data.begin();
	for(; iter != data.end(); iter++)
	{
		if(iter->first <= max)
			fprintf(fp, "%d\t%d\n", iter->first, iter->second);
		else
			maxCount += iter->second;
	}
	fprintf(fp, ">%d\t%d\n", max, maxCount);

}

