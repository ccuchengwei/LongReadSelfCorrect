#ifndef Kmer_H
#define Kmer_H

#include "BWTIndexSet.h"
#include "BWTAlgorithms.h"
#include "Alphabet.h"
#include <algorithm>
#include <map>
#include <memory>
//Implemntaiont of kmer object on each position of some sequence;
//theoretically all kmers could be established on previous smaller one. Noted by KuanWeiLee 18/3/12
/*
ex: On position 14, 4-mer is established first, and 6-mer is established on 4-mer, so does 10-mer on 4-mer.
              14(pos)
			  .
===========================================================================> sequence
              ---> size == 4
			  -----> size == 6 
			  ---------> size == 10
*/
class KmerFeature
{
	public:

		static thread_local std::map<int, std::unique_ptr<KmerFeature[]> > kmerRec;
	
		//Null kmer
		KmerFeature()
		:	count(nullptr),
			indices(nullptr){ }
	
		//Copy kmer
		KmerFeature(const KmerFeature& base)
		:	count(base.getCount()),
			indices(base.getIndex()),
			word(base.getWord()),
			size(base.getSize()),
			biInterval(base.getInterval()),
			isFake(base.getPseudo()),
			frequency(base.getFreq()){ }
		
		//Base(or not) kmer
		KmerFeature(
				const BWTIndexSet* indexSet,
				const std::string & seq,
				size_t pos,
				int len,
				const KmerFeature* base = nullptr)
		:	count(base == nullptr ? std::unique_ptr<int[]>(new int[DNA_ALPHABET::size]{0}) : (*base).getCount()),
			indices(base == nullptr ? indexSet : (*base).getIndex()),
			word(base == nullptr ?  seq.substr(pos, len) : (*base).getWord()),
			size(base == nullptr ? word.length() : (*base).getSize()),
			biInterval(base == nullptr ? BWTAlgorithms::findBiInterval(*indices, word, count.get()) : (*base).getInterval())
		{
			size_t seqlen = seq.length();
			for(size_t i = (pos + size); (i < seqlen) && (size < len); i++)
			{
				char b = seq[i];
				expand(b);
			}
			isFake = (len != size);
			frequency = biInterval.getFreq();
		}
	
		~KmerFeature(){ }
	
		inline KmerFeature& operator=(const KmerFeature& other)
		{
		//	if(&other == this) return *this;
			count = other.getCount();
			indices = other.getIndex();
			word = other.getWord();
			size = other.getSize();
			biInterval = other.getInterval();
			isFake = other.getPseudo();
			frequency = other.getFreq();
			return *this;
		}
	
		inline std::unique_ptr<int[]> getCount() const
		{
			const int* oldOne = count.get();
			std::unique_ptr<int[]> newOne(new int[DNA_ALPHABET::size]);
			std::copy(oldOne, (oldOne + 4), newOne.get());
			return newOne;
		}
		inline const BWTIndexSet* getIndex() const { return indices; }
		inline std::string getWord() const { return word; }
		inline int getSize() const { return size; }
		inline BiBWTInterval getInterval() const { return biInterval; }
		inline bool getPseudo() const { return isFake; }
		inline int getFreq() const { return isFake ? 0 : frequency; }
	
		inline bool isValid() const { return biInterval.isValid(); }
		inline void expand(char b)
		{
			if(b == 0) return;
			size++;
			word += b;
			BWTAlgorithms::updateBiInterval(biInterval, b, *indices, count.get());
			frequency = biInterval.getFreq();
		}
		inline void shrink(int len, bool update = false)
		{
			size -= len;
			for(std::string::iterator iter = (word.begin() + size); iter != word.end(); iter++)
				count[DNA_ALPHABET::getIdx(*iter)]--;
			word.erase(size, len);
			if(!update) return;
			biInterval = BWTAlgorithms::findBiInterval(*indices, word);
			frequency = biInterval.getFreq();
		}
		inline bool isLowComplexity(float t = 0.7) const
		{
			bool isMonmer = false, isDimer = false;
			int num = 0;
			for(int i = 0; i < DNA_ALPHABET::size; i++)
			{
				isMonmer = isMonmer || ((float)count[i]/size >= t);
				num += (count[i] == 0 ? 1 : 0);
			}
			isDimer = (num == 2);
			return isMonmer || isDimer;
		}

	private:
		std::unique_ptr<int[]> count;
		const BWTIndexSet* indices;
		std::string word;
		int size;
		BiBWTInterval biInterval;
		//'isFake' is set only when initialized.
		bool isFake;
		int frequency;
	
};

#endif