// CheckKmerProcess - Check kmer distribution and error condition
//

#ifndef CHECKKMERPROCESS_H
#define CHECKKMERPROCESS_H

#include <iostream>
#include <map>
#include <utility>
#include <cstdio>
#include "BWTIndexSet.h"
#include "SequenceWorkItem.h"
#include "KmerDistribution.h"

class CodeBlock
{
	public:
		CodeBlock(
				int _start,
				int _end,
				const std::string& _code,
				bool _rvc)
		:	start(_start),
			end  (_end),
			code (_code),
			rvc  (_rvc){ }
		~CodeBlock() = default;
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

//Helpful typedef
typedef std::map<int, std::ostream*> ostreamPtrMap;
typedef std::map<int, KmerDistribution> kdMap;
// Parameter object
struct CheckKmerParameters
{
	BWTIndexSet indices;
	std::string directory;
	std::pair<int ,int> size;
	std::map<std::string, std::list<CodeBlock> >* pAlignLog;
	//mode(true/false) ? seq : seed
	bool mode;
};

struct CheckKmerResult
{
	kdMap corKdMap;
	kdMap errKdMap;
	//0:none 1:true 2:false
	std::string readid;
	int status[3]{0, 0 ,0};
};

//
class CheckKmerProcess
{
	public:
		CheckKmerProcess(CheckKmerParameters params):m_params(params){ }
		~CheckKmerProcess(){ }
		CheckKmerResult process(const SequenceWorkItem& workItem);
	
	private:
		void scan(int ksize, const CodeBlock& block, const std::string& query, CheckKmerResult& result);
		std::string fetch(const std::string& in, int pos, int step);
		int sum(const std::string& in);
		int getPys(int pos, int len);
		bool validate(int pos, int ksize, const CodeBlock& block, const std::string& query);
		CheckKmerParameters m_params;
};
//
class CheckKmerPostProcess
{
	public:
		CheckKmerPostProcess(CheckKmerParameters params);
		~CheckKmerPostProcess();

		void process(const SequenceWorkItem& workItem, const CheckKmerResult& result);

	private:
		CheckKmerParameters m_params;
		kdMap m_corKdMap;
		kdMap m_errKdMap;
		ostreamPtrMap m_pCorWriterMap;
		ostreamPtrMap m_pErrWriterMap;
		int* m_status;
		FILE *m_pStatusWriter;
};

#endif
