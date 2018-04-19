#ifndef __INTERVAL_TREE_H
#define __INTERVAL_TREE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>

template <class T, typename K = std::size_t>
class TreeInterval
{
	public:
	    TreeInterval(
				K s,
				K e,
				const T& v)
		:	start(s),
			stop(e),
			value(v){ }
		K start;
		K stop;
		T value;
};

template <class T, typename K>
K intervalStart(const TreeInterval<T,K>& i)
{
    return i.start;
}

template <class T, typename K>
K intervalStop(const TreeInterval<T,K>& i)
{
    return i.stop;
}

template <class T, typename K>
std::ostream& operator<<(std::ostream& out, TreeInterval<T,K>& i)
{
    out << "TreeInterval(" << i.start << ", " << i.stop << "): " << i.value;
    return out;
}

template <class T, typename K = std::size_t>
class IntervalStartSorter
{
	public:
		bool operator() (const TreeInterval<T,K>& a, const TreeInterval<T,K>& b)
		{
			return a.start < b.start;
		}
};

template <class T, typename K = std::size_t>
class IntervalTree
{
	public:
		typedef TreeInterval<T,K> interval;
		typedef std::vector<interval> intervalVector;
		typedef IntervalTree<T,K> intervalTree;
	
		IntervalTree(void)
		:	left(nullptr),
			right(nullptr),
			center(0){ }
	
		IntervalTree(const intervalTree& other)
		:   intervals(other.intervals),
			left(other.left ? copyTree(*other.left) : nullptr),
			right(other.right ? copyTree(*other.right) : nullptr),
			center(other.center){ }
	
		IntervalTree(
				intervalVector& ivals,
				std::size_t depth = 16,
				std::size_t minbucket = 64,
				K leftextent = 0,
				K rightextent = 0,
				std::size_t maxbucket = 512);
	
		~IntervalTree(void) = default;
	
		IntervalTree<T,K>& operator=(const intervalTree& other);
		intervalVector findOverlapping(K start, K stop) const;
		void findOverlapping(K start, K stop, intervalVector& overlapping) const;
		intervalVector findContained(K start, K stop) const;
		void findContained(K start, K stop, intervalVector& contained) const;	
	
		intervalVector intervals;
		std::unique_ptr<intervalTree> left;
		std::unique_ptr<intervalTree> right;
		K center;

	private:
		std::unique_ptr<intervalTree> copyTree(const intervalTree& orig);
};

#endif
