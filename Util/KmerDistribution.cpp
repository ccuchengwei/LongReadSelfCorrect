//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// KmerDistribution - Histogram of kmer frequencies
//
#include <algorithm>
#include <cmath>
#include <iostream>
#include "KmerDistribution.h"

KmerDistribution::KmerDistribution()
:	total(0),
	q1(0),
	q2(0),
	q3(0),
	min(0),
	max(0),
	mode(0),
	sdv(0),
	repeatKmerCutoff(0){ }

void KmerDistribution::operator+=(const KmerDistribution& other)
{
	for(const auto& iter : other.data)
		data[iter.first] += iter.second;
	total += other.total;
}

int KmerDistribution::getNumberWithCount(int n) const
{
	std::map<int, int>::const_iterator iter = data.find(n);
	return iter != data.end() ? iter->second : 0;
}

int KmerDistribution::getQuartile(int n) const
{
	switch(n)
	{
		case 1: return q1;
		case 2: return q2;
		case 3: return q3;
		default:
			std::cerr << "Quartile option: 1,2,3\n";
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

void KmerDistribution::computeKDAttributes()
{
	int low = total*1/4;
	int mid = total*2/4;
	int upp = total*3/4;
	int prev = 0;
	int curr = 0;
	int most = 0;
	for(const auto& iter : data)
	{
		if(iter.second > most)
		{
			most = iter.second;
			mode = iter.first;
		}
		
		prev = curr;
		curr += iter.second;
		if(low >= prev && low <= curr) q1 = iter.first;
		if(mid >= prev && mid <= curr) q2 = iter.first;
		if(upp >= prev && upp <= curr) q3 = iter.first;
	//	if(q3 > 0) break;
	}
	
	int iqr = q3 - q1;
	int small = q1 - (int)(iqr*1.5);
	int large = q3 + (int)(iqr*1.5);
	prev = curr = 0;
	for(const auto& iter :data)
	{
		prev = curr;
		curr = iter.first;
		if(min == 0 && curr >= small) min = curr;
		if(prev <= large && curr > large) max = prev;
	//	if(max > 0) break;
	}
	if(max == 0) max = curr;
	
	int sqsum = 0;
	std::for_each(data.begin(), data.end(), [&](std::pair<int,int> x)mutable{sqsum += x.second*pow((x.first - q2), 2);});
	double variance = (double)sqsum/(total - 1);
	sdv = sqrt(variance);
	
	// double freq95 = getCutoffForProportion(0.95);
	// repeatKmerCutoff = sdv > q2*2? q2*1.5: q2*1.3;
	repeatKmerCutoff = q2*1.3;
	// repeatKmerCutoff = getCutoffForProportion(0.8);
	// repeatKmerCutoff = (double) q2*(0.39+0.53* (freq95/(double)q2));
}

std::ostream& operator<<(std::ostream& out, const KmerDistribution& o)
{
	out << o.min << ' ' << o.q1 << ' ' << o.q2 << ' ' << o.q3 << ' ' << o.max;
	return out;
}

void compare(std::ostream& t, std::ostream& v, int cov, int ksize, KmerDistribution& c, KmerDistribution& e)
{
	c.computeKDAttributes();
	e.computeKDAttributes();
	
	t << cov << ' ' << ksize << " | " << e << " | " << c << '\n';
	
	int value = 0;
	if(c.min >= e.max) value = c.min;
	else if(c.q1 >= e.q3) value = c.q1;
	else value = c.q1;
	v << cov << ' ' << ksize << ' ' << value << '\n';
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
/***********/