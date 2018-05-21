#ifndef BARCODE_H
#define BARCODE_H

#include <map>
#include <vector>

extern const std::map<char, int> base_hex;
extern const std::map<char, int> char_int;

class BCode
{
	public:
		BCode(
				int s,
				int e,
				const std::string& c,
				bool r)
		:	start(s),
			end  (e),
			code (c),
			rvc  (r) { }
		~BCode() = default;
		
		typedef std::vector<BCode> BCodeVector;
		static std::map<std::string, BCodeVector>& Log();
		static void load(const std::string& barcodefile);
		static std::string fetch(const std::string& in, int pos, int step);
		static int sum(const std::string& in);
		static int getPys(int pos, int len);
		static bool validate(int pos, int ksize, const BCode& block, const std::string& seq);
		
		inline int getStart() const { return start; }
		inline int getEnd  () const { return end;   }
		inline const std::string& getCode() const { return code; }
		inline bool getRvc() const { return rvc; }
		
	private:
		int start;
		int end;
		std::string code;
		bool rvc;
};

#endif