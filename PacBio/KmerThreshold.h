#ifndef KMERTHRESHOLDTABLE_H
#define KMERTHRESHOLDTABLE_H

class KmerThreshold
{
	public:
		KmerThreshold(const KmerThreshold&) = delete;
		void operator=(const KmerThreshold&) = delete;
		
		void initialize(int s, int e, int c, const std::string& d);
		void write(std::ostream& out = std::cout);
		
		inline static KmerThreshold& Instance()
		{
			static KmerThreshold instance;
			return instance;
		}
		inline float get(int mode, int ksize)
		{
			assert(mode >= 0 && mode <= 2);
			assert(ksize >= 0 && ksize <= (end + 1));
			return table[mode][ksize];
		}
		
	private:
		KmerThreshold();
		~KmerThreshold();
		
		float calculate(int mode, int x, int y);
		
		int start;
		int end;
		int cov;
		float* table[3];
		std::ostream* pTableWriter;
		bool set;
};

#endif