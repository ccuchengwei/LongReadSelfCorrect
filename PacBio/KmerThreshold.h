#ifndef KMERTHRESHOLDTABLE_H
#define KMERTHRESHOLDTABLE_H

class KmerThreshold
{
	public:
		KmerThreshold(const KmerThreshold&) = delete;
		void operator=(const KmerThreshold&) = delete;
		
		void set(int _start, int _end, int _cov, const std::string& dir);
		void write(std::ostream& outfile = std::cout);
		
		inline static KmerThreshold& Instance()
		{
			static KmerThreshold instance;
			return instance;
		}
		inline float get(int type, int ksize)
		{
			assert(type >= 0 && type <= 2);
			assert(ksize >= 0 && ksize <= (end + 1));
			return table[type][ksize];
		}
		
	private:
		KmerThreshold();
		~KmerThreshold();
		
		float calculate(int type, int x, int y);
		
		int start;
		int end;
		int cov;
		float* table[3];
		std::ostream* pTableWriter;
};

#endif