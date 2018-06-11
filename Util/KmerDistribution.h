//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// KmerDistribution - Histogram of kmer frequencies
//
#ifndef KMERDISTRIBUTION_H
#define KMERDISTRIBUTION_H

#include <vector>
#include <map>
#include <cstdio>
#include <iostream>

class KmerDistribution
{
	public:
		
		KmerDistribution();
		~KmerDistribution(void) = default;
		void operator+=(const KmerDistribution& other);
		
		inline void add(int kmerFreqs){ data[kmerFreqs]++; total++; }
		inline double getSdv() const{ return sdv; }
		inline int getTotalKmers() const{ return total; }
		inline int getMode() const{ return mode; }
		inline int getRepeatKmerCutoff() const{ return repeatKmerCutoff; }
		
		int getNumberWithCount (int n) const;
		int getQuartile(int n) const;
		
		// Returns the proportion of the distribution less than or equal to n
		double getCumulativeProportionLEQ(int n) const;
		
		// Returns the smallest value n such that the proportion of the data less
		// than n is p. This is the inverse of getCumulativeProportionLEQ.
		int getCutoffForProportion(double p) const;
		
		//compute quartiles ,mode and sdv
		void computeKDAttributes();
		
		//Write data to file
		friend std::ostream& operator<<(std::ostream& out, const KmerDistribution& o);
		
		friend void compare(std::ostream& t, std::ostream& v, int cov, int ksize, KmerDistribution& c, KmerDistribution& e);
		
		//Legacy Part
		/***********/
		// Returns the predicted boundary of the left tail which holds erroneous kmers
		int findFirstLocalMinimum() const;
		int findErrorBoundary() const;
		int findErrorBoundaryByRatio(double ratio) const;
		
		// Get the mode of the distribution if the first n values are ignored
		int getCensoredMode(int n) const;
		std::vector<int> toCountVector(int max) const;
		void print(int max) const; 
		void print(FILE *file, int max) const;
		/***********/
		
	private:
		std::map<int, int> data;
		int total;
		int q1;
		int q2;
		int q3;
		int min;
		int max;
		int mode;
		double sdv;
		int repeatKmerCutoff;
};

#endif
