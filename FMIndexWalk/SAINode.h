//----------------------------------------------
// Copyright 2015 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------

// 
// The search tree represents a traversal through implicit FM-index graph
//
#ifndef SAINODE_H
#define SAINODE_H

#include <list>
#include "BWT.h"
#include "BWTAlgorithms.h"
class SagNode
{
    public: 
    SagNode(BWTInterval interval,size_t index,std::string kmerstring):m_interval(interval),m_index(index),kmerstring(kmerstring){};
    SagNode(BWTInterval interval,size_t index,std::string kmerstring,SagNode* parentnode,size_t value):m_interval(interval),m_index(index),kmerstring(kmerstring)

    {
        makeinedge(parentnode,value);
        parentnode->makeoutedge(this,value);
    };
    ~SagNode();
    void makeoutedge(SagNode* childnode,size_t value){ outedge.push_back(std::make_pair(childnode,value));};
    void makeinedge(SagNode* parentnode,size_t value){ inedge.push_back(std::make_pair(parentnode,value)); }; 
            
            
     
    BWTInterval m_interval;
    size_t m_index;
    std::string kmerstring; 
    std::vector<std::pair<SagNode*,size_t>> outedge,inedge;
    private:
    
   
    
};
class SAINode;
typedef std::list<SAINode*> SAINodePtrList;

// Base class represent search nodes for all sorts of applications using FM-index walk
class SAINode
{
    public:

        //
        // Functions
        //
        SAINode(const std::string* pQuery, SAINode* parent);
        ~SAINode();

        // Add a child node to this node with the given label
        // Returns a pointer to the created node
        SAINode* createChild(const std::string& label);

        // Extend the label of this node by l
        void extend(const std::string& ext);

        // Return a suffix of length l of the string represented by this node
        std::string getSuffix(size_t l) const;

        // Return the complete sequence of the string represented by the branch
        std::string getFullString() const;

        // Initialize or update the alignment data
        //void computeInitialAlignment(const std::string& initialLabel, int queryAlignmentEnd, int bandwidth);
        void computeInitial(const std::string& initialLabel);

        // Recursive function to print all the strings represented
        // by this node and all its children.
        void printAllStrings(const std::string& parent) const;

        size_t getKmerCount(){return m_totalKmerCount;};
		
		size_t getLastKmerCount(){return m_lastKmerCount;};
		
        void addKmerCount(size_t currKmerCount){ 
			m_totalKmerCount += currKmerCount; 
			m_lastKmerCount = currKmerCount;
		};

        void reduceKmerCount(size_t currKmerCount){ 
			m_totalKmerCount -= currKmerCount; 
			// m_lastKmerCount = currKmerCount;
		};
		
		// get number of children
		size_t getNumOfChild(){return m_children.size();}
		
		// kill all children
		void killAllChild()
		{
			m_children.clear();
		}
		
		// trim label
		void trimLabel(size_t trimLen)
		{
			m_label = m_label.substr(0, m_label.length()-trimLen);
		}
		
		inline void setUpdate(bool flag)
		{
			m_updated = flag;
		}
		
		inline bool isUpdated(){return m_updated;};


    private:

        //
        // Data
        //

        // The extension string from the parent
        std::string m_label;
       size_t m_totalKmerCount;
		size_t m_lastKmerCount;
		bool m_updated;
        // size_t m_kmersize;


	protected:
        // The query string being threaded through the graph
        const std::string* m_pQuery;

        // The parent node, can be NULL
        SAINode* m_pParent;
        SAINodePtrList m_children;

};


//
// SAIntervalNode for implementation of FM-index walk between paired end reads or validation
//
class SAIntervalNode : public SAINode
{
    public:

        //
        // Functions
        //
        SAIntervalNode(const std::string* pQuery, SAIntervalNode* parent):SAINode(pQuery,parent)
		{
		}
        ~SAIntervalNode(){};

		// Add a child node to this node with the given label
        // Returns a pointer to the created node
        SAIntervalNode* createChild(const std::string& label);

        BWTInterval fwdInterval;
        BWTInterval rvcInterval;

};

// leaves of SAIntervalNode
typedef std::list<SAIntervalNode*> STNodePtrList;

// Object to hold the result of the extension process
struct SAIntervalNodeResult
{
    std::string thread;
	size_t SAICoverage;
	int SAIntervalSize;
    double errorRate;
};
typedef std::vector<SAIntervalNodeResult> SAIntervalNodeResultVector;



//
// SAIOverlapNode for implementation of overlap computation using FM-index walk
//
class SAIOverlapNode : public SAINode
{
    public:

        //
        // Functions
        //
        SAIOverlapNode(const std::string* pQuery, SAIOverlapNode* parent):SAINode(pQuery,parent)
		{
			lastSeedIdx=totalSeeds=lastOverlapLen=currOverlapLen=queryOverlapLen=numOfErrors=0;
			lastSeedIdxOffset=0;
		}
        ~SAIOverlapNode(){};

		// Add a child node to this node with the given label
        // Returns a pointer to the created node
        SAIOverlapNode* createChild(const std::string& label);

        BWTIntervalPair currIntervalPair;
		
		// last matched seed index
		size_t lastSeedIdx;
		// last overlap length when matching last seed
		size_t lastOverlapLen;
		size_t totalSeeds;
		// current overlap length on the subject increases wrt each FM-index extension
		size_t currOverlapLen;
		// number of SNPs or indels
		size_t numOfErrors;
		// index offset to the center
		int lastSeedIdxOffset;

		// index of the init seed
		int initSeedIdx;
		
		// current overlap length on the query
		size_t queryOverlapLen;
		
		
		SagNode* currseedindex;
};

// leaves of SAIOverlapNode
typedef std::list<SAIOverlapNode*> SONodePtrList;

//
// SAIOverlapNode2 for implementation of overlap computation using FM-index walk
// It's used by ShortReadOverlapTree, v20160802 by Ya.
//
class SAIOverlapNode2 : public SAINode
{
   public:
        //
        // Functions
        //
        SAIOverlapNode2(const std::string* pQuery, SAIOverlapNode2* parent):SAINode(pQuery,parent)
		{
			lastSeedIdx=totalSeeds=lastOverlapLen=currOverlapLen=queryOverlapLen=numOfErrors=0;
			numRedeemSeed=0;
			lastSeedIdxOffset=0;
		}
        ~SAIOverlapNode2(){};

		// Add a child node to this node with the given label
        // Returns a pointer to the created node
        SAIOverlapNode2* createChild(const std::string& label);

        BWTInterval fwdInterval;
        BWTInterval rvcInterval;
		
		// last matched seed index
		size_t lastSeedIdx;
		// error seed begin idx
		// size_t errorSeedBeginIdx;
		// number of redeem seeds
		double numRedeemSeed;
		// last overlap length when matching last seed
		size_t lastOverlapLen;
		size_t totalSeeds;
		// current overlap length on the subject increases wrt each FM-index extension
		size_t currOverlapLen;
		// number of SNPs or indels
		size_t numOfErrors;
		// index offset to the center
		int lastSeedIdxOffset;
		// index of the init seed
		int initSeedIdx;
		// current overlap length on the query
		size_t queryOverlapLen;


        

};

// leaves of SAIOverlapNode2
typedef std::list<SAIOverlapNode2*> SONode2PtrList;

class SAIOverlapNode3;
// leaves of SAIOverlapNode3
typedef std::list<SAIOverlapNode3*> SONode3PtrList;


//
// SAIOverlapNode2 for implementation of overlap computation using FM-index walk
// It's used by ShortReadOverlapTree, v20160802 by Ya.
//
class SAIOverlapNode3 : public SAINode
{
    public:
        //
        // Functions
        //
        SAIOverlapNode3(const std::string* pQuery, SAIOverlapNode3* parent):SAINode(pQuery,parent)
		{
			lastSeedIdx=totalSeeds=lastOverlapLen=currOverlapLen=queryOverlapLen=numOfErrors=0;
			numRedeemSeed=0;
			lastSeedIdxOffset=0;
		}
        ~SAIOverlapNode3();

		// Add a child node to this node with the given label
        // Returns a pointer to the created node
        SAIOverlapNode3* createChild(const std::string& label);


        BWTInterval fwdInterval;
        BWTInterval rvcInterval;
		
		// last matched seed index
		size_t lastSeedIdx;
		// error seed begin idx
		// size_t errorSeedBeginIdx;
		// number of redeem seeds
		double numRedeemSeed;
		// last overlap length when matching last seed
		size_t lastOverlapLen;
		size_t totalSeeds;
		// current overlap length on the subject increases wrt each FM-index extension
		size_t currOverlapLen;
		// number of SNPs or indels
		size_t numOfErrors;
		// index offset to the center
		int lastSeedIdxOffset;
		// index of the init seed
		int initSeedIdx;
		// current overlap length on the query
		size_t queryOverlapLen;
        // index of the result and index of the matchpoint
        std::pair <int,int> resultindex = std::make_pair(-1,-1);
        // size_t currkmersize;
        
        int WeightofContinuousChar;
        std::vector<double> LocalErrorRateRecord;
        std::vector<double> GlobalErrorRateRecord;
        SONode3PtrList m_children3;
  
        
        
        
};

// leaves of SAIOverlapNode3
typedef std::list<SAIOverlapNode3*> SONode3PtrList;


#endif
