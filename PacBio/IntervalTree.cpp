#include "IntervalTree.h"
#include "IntervalTreeInstantiation.cpp"

template <class T, typename K>
IntervalTree<T,K>::IntervalTree(
		intervalVector& ivals,
		std::size_t depth,
		std::size_t minbucket,
		K leftextent,
		K rightextent,
		std::size_t maxbucket)
:	left(nullptr),
	right(nullptr),
	center(0)
{
	K leftp = leftextent, rightp = rightextent, centerp = 0;
	if(leftp == 0 && rightp == 0)
		std::sort(ivals.begin(), ivals.end(), std::greater<interval>());
	
	if(--depth == 0 || ivals.size() < minbucket)
		intervals = ivals;
	else
	{
		leftp = ivals.back().start;
		rightp = max_element(ivals.begin(), ivals.end())->stop;
		
		centerp = ivals[ivals.size() >> 1].start;
		center = centerp;

		intervalVector lefts;
		intervalVector rights;

		for(const auto& interval : ivals)
		{
			if(interval.stop < center)
				lefts.push_back(interval);
			else if(interval.start > center)
				rights.push_back(interval);
			else
				intervals.push_back(interval);
		}

		if (!lefts.empty())
			left = std::unique_ptr<intervalTree>(new intervalTree(lefts, depth, minbucket, leftp, centerp));
		if (!rights.empty())
			right = std::unique_ptr<intervalTree>(new intervalTree(rights, depth, minbucket, centerp, rightp));
	}
}

template <class T, typename K>
typename IntervalTree<T,K>::intervalTree& IntervalTree<T,K>::operator=
(const typename IntervalTree<T,K>::intervalTree& other)
{
	center = other.center;
	intervals = other.intervals;
	left = other.left ? copyTree(*other.left) : nullptr;
	right = other.right ? copyTree(*other.right) : nullptr;
	return *this;
}
//Theorectically speaking, one bi(fwd/rvc) interval would (overlap/be contained in) another,
//if these two strings have common (prefix/suffix); it's impossible to find two realistic intervals
//'interleave' each other. Noted by KuanWeiLee 18/4/10
template <class T, typename K>
typename IntervalTree<T,K>::intervalVector IntervalTree<T,K>::findOverlapping
(K start, K stop) const
{
	intervalVector ov;
	this->findOverlapping(start, stop, ov);
	return ov;
}

template <class T, typename K>
void IntervalTree<T,K>::findOverlapping
(K start, K stop, intervalVector& overlapping) const
{
	if(!intervals.empty() && ! (stop < intervals.back().start))
	{
		for(const auto& interval : intervals)
		{
			if (interval.start <= start && interval.stop >= stop)
				overlapping.push_back(interval);
		}
	}

	if(left && start < center)
		left->findOverlapping(start, stop, overlapping);

	if(right && stop > center)
		right->findOverlapping(start, stop, overlapping);

}

template <class T, typename K>
typename IntervalTree<T,K>::intervalVector IntervalTree<T,K>::findContained
(K start, K stop) const
{
	intervalVector contained;
	this->findContained(start, stop, contained);
	return contained;
}

template <class T, typename K>
void IntervalTree<T,K>::findContained
(K start, K stop, intervalVector& contained) const
{
	if(!intervals.empty() && ! (stop < intervals.back().start))
	{
		for(const auto& interval : intervals)
		{
			if(interval.start >= start && interval.stop <= stop)
				contained.push_back(interval);
		}
	}

	if(left && start <= center)
		left->findContained(start, stop, contained);

	if(right && stop >= center)
		right->findContained(start, stop, contained);

}

template <class T, typename K>
std::unique_ptr<typename IntervalTree<T,K>::intervalTree> IntervalTree<T,K>::copyTree
(const typename IntervalTree<T,K>::intervalTree& orig)
{
	return std::unique_ptr<intervalTree>(new intervalTree(orig));
}
