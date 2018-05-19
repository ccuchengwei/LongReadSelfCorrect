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
#include <iostream>

class KmerDistribution
{
	public:
		
		KmerDistribution(void) = default;
		~KmerDistribution(void) = default;
		void operator+=(const KmerDistribution& temp);
		
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
		void computeKDAttributes(float censor = 0.0f);
		
		//Write data to file
		void write(std::ostream& out, int mode = 0) const;
		
		
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
		std::string readid = "";
		int total = 0;
		int median = 0;
		int first_quartile = 0;
		int third_quartile = 0;
		int mode = 0;
		double sdv = 0;
		int repeatKmerCutoff = 0;
};

#endif
