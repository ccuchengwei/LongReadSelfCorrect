#ifndef Kmer_H
#define Kmer_H

#include <algorithm>
#include <map>
#include <memory>
#include "BWTIndexSet.h"
#include "BWTAlgorithms.h"
#include "Alphabet.h"

/*
Implemntaiont of kmer object on each position of some sequence;
theoretically all kmers could be established on previous smaller one. Noted by KuanWeiLee 18/3/12
ex: On position 14, 4-mer is established first, and 6-mer is established on 4-mer, so does 10-mer on 6-mer.
			. 14[pos]
===========================================================================> sequence
            ---> size == 4
			-----> size == 6 
			---------> size == 10
*/
class KmerFeature
{
	public:
		inline static std::map<int ,std::unique_ptr<KmerFeature[]> >& Log()
		{
			static thread_local std::map<int, std::unique_ptr<KmerFeature[]> > log;
			return log;
		}
	
		KmerFeature(void) = default;
		~KmerFeature(void) = default;
		
		//Copy kmer
		KmerFeature(const KmerFeature& base){ *this = base; }
		
		//Base(or not) kmer
		KmerFeature(
				const BWTIndexSet& indices,
				const std::string& seq,
				size_t pos,
				int len,
				const KmerFeature* base = nullptr)
		{
			if(base == nullptr)
			{
				this->count      = std::unique_ptr<int[]>(new int[DNA_ALPHABET::size]{0});
				this->indices    = indices;
				this->word       = seq.substr(pos, len);
				this->size       = word.length();
				this->biInterval = BWTAlgorithms::findBiInterval(this->indices, this->word, this->count.get());
			}
			else
			{
				*this = *base;
				assert(this->size < len);
				for(size_t i = (pos + this->size); (i < seq.length()) && (this->size < len); i++)
				{
					char b = seq[i];
					expand(b);
				}
			}
			this->fake = (len != this->size);
			this->frequency = this->biInterval.getFreq();
		}
	
		inline KmerFeature& operator=(const KmerFeature& other)
		{
			if(&other == this) return *this;
			
			this->count      = other.getCount();
			this->indices    = other.indices;
			this->word       = other.word;
			this->size       = other.size;
			this->biInterval = other.biInterval;
			this->fake       = other.fake;
			this->frequency  = other.frequency;
			return *this;
		}
	
		inline std::unique_ptr<int[]> getCount() const
		{
			const int* orig = this->count.get();
			std::unique_ptr<int[]> copy(new int[DNA_ALPHABET::size]);
			std::copy(orig, (orig + 4), copy.get());
			return copy;
		}
		
		inline std::string getWord() const { return this->word; }
		inline int getSize() const { return this->size; }
		inline int getFreq() const { return this->fake ? -1 : this->frequency; }
	
		inline void expand(char b)
		{
			if(b == 0) return;
			this->size++;
			this->word += b;
			BWTAlgorithms::updateBiInterval(this->biInterval, b, this->indices, this->count.get());
			this->frequency = this->biInterval.getFreq();
		}
		
		inline void shrink(int len, bool update = false)
		{
			assert(len < this->size);
			this->size -= len;
			for(std::string::iterator iter = (this->word.begin() + this->size); iter != this->word.end(); iter++)
				this->count[DNA_ALPHABET::getIdx(*iter)]--;
			this->word.erase(this->size, len);
			if(!update) return;
			this->biInterval = BWTAlgorithms::findBiInterval(this->indices, this->word);
			this->frequency = this->biInterval.getFreq();
		}
		
		inline bool isFake() const { return this->fake; };
		inline bool isValid() const { return this->biInterval.isValid(); }
		
		inline bool isLowComplexity(float m = 0.7, float d = 0.9) const
		{
			const int* orig = this->count.get();
			int* copy = new int[4];
			std::copy(orig, (orig + 4), copy);
			std::sort(copy, (copy + 4));
			bool isMonmer = (float)copy[3]/this->size >= m;
			bool isDimer  = (float)(copy[2] + copy[3])/this->size >= d;
			delete[] copy;
			return isMonmer || isDimer;
		}

	private:
		std::unique_ptr<int[]> count;
		BWTIndexSet indices;
		std::string word;
		int size;
		BiBWTInterval biInterval;
		bool fake; //'fake' is set only when initialized.
		int frequency;
};

#endif