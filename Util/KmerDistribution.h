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
#include <cstddef>
#include <iostream>

class KmerDistribution
{
    public:
		
        KmerDistribution(std::string name = ""):m_readid(name){};
        ~KmerDistribution(){};
        inline void add(int kmerFreqs){m_data[kmerFreqs]++; m_total++;};
		inline size_t getQuartile(int n) const
		{
			switch(n)
			{
				case 1: return m_first_quartile;
				case 2: return m_median;
				case 3: return m_third_quartile;
				default: 
						std::cout << "Quartile options : 1 , 2 ,3\n";
						exit(EXIT_FAILURE);
			}
		};
		inline double getSdv() const{return m_std;};
		inline size_t getTotalKmers() const{return m_total;};
		inline size_t getMode() const{return m_mode;};
		inline size_t getRepeatKmerCutoff() const{return m_repeatKmerCutoff;};
		inline size_t getNumberWithCount (size_t n) const
		{
			iteratorKmerFreqsMap iterElement = m_data.find(n);
			return (iterElement != m_data.end() ? iterElement->second : 0);
		};				
        // Returns the proportion of the distribution less than or equal to n
        double getCumulativeProportionLEQ(int n) const;
        // Returns the smallest value n such that the proportion of the data less
        // than n is p. This is the inverse of getCumulativeProportionLEQ.
        size_t getCutoffForProportion(double p) const;
		//compute quartiles ,mode and stdev
		void computeKDAttributes(float censor = 0.0f);
		//Write data to file
        enum TYPE {DATA,ATTRIBUTE};
		void write(std::ostream& outfile,TYPE mode = DATA) const;
		
		
		//Below are legacy codes. Noted by KuanWeiLee 20171027
		/***********************************************************************************/
		// Returns the predicted boundary of the left tail which holds erroneous kmers
		int findFirstLocalMinimum() const;
        int findErrorBoundary() const;
        int findErrorBoundaryByRatio(double ratio) const;
		// Get the mode of the distribution if the first n values are ignored
        size_t getCensoredMode(size_t n) const;
        std::vector<int> toCountVector(int max) const;
        void print(int max) const; 
        void print(FILE* file, int max) const; 
    private:

		//kmerFreqs --> number
		typedef std::map<size_t,size_t> KmerFreqsMap;
		typedef KmerFreqsMap::const_iterator iteratorKmerFreqsMap;
		KmerFreqsMap m_data;
		std::string m_readid;
		size_t m_total = 0;
		size_t m_median = 0;
		size_t m_first_quartile = 0;
		size_t m_third_quartile = 0;
		size_t m_mode = 0;
		double m_std = 0;
		size_t m_repeatKmerCutoff = 0;
};

#endif
