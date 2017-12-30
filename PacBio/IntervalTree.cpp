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
	{
		IntervalStartSorter<T,K> intervalStartSorter;
		std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
	}
	
	if(--depth == 0 || ivals.size() < minbucket)
		intervals = ivals;
	else
	{
		leftp = ivals.front().start;
		std::vector<K> stops;
		stops.resize(ivals.size());
		transform(ivals.begin(), ivals.end(), stops.begin(), intervalStop<T,K>);
		rightp = *max_element(stops.begin(), stops.end());
		
		//centerp = (leftp + rightp) >> 1;
		centerp = ivals[ivals.size() >> 1].start;
		center = centerp;

		intervalVector lefts;
		intervalVector rights;

		for (typename intervalVector::const_iterator i = ivals.begin(); i != ivals.end(); ++i)
		{
			const interval& interval = *i;
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
	
	/*
	--depth;
	IntervalStartSorter<T,K> intervalStartSorter;
	if (depth == 0 || ivals.size() < minbucket )
	{
		std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
		intervals = ivals;
	} 
	else
	{
		if (leftextent == 0 && rightextent == 0)
		{
			// sort intervals by start
			std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
		}

		K leftp = 0;
		K rightp = 0;
		K centerp = 0;

		if (leftextent || rightextent)
		{
			leftp = leftextent;
			rightp = rightextent;
		}
		else
		{
			leftp = ivals.front().start;
			std::vector<K> stops;
			stops.resize(ivals.size());
			transform(ivals.begin(), ivals.end(), stops.begin(), intervalStop<T,K>);
			rightp = *max_element(stops.begin(), stops.end());
		}

		//centerp = ( leftp + rightp ) / 2;
		centerp = ivals.at(ivals.size() / 2).start;
		center = centerp;

		intervalVector lefts;
		intervalVector rights;

		for (typename intervalVector::const_iterator i = ivals.begin(); i != ivals.end(); ++i)
		{
			const interval& interval = *i;
			if(interval.stop < center)
			{
				lefts.push_back(interval);
			}
			else if(interval.start > center)
			{
				rights.push_back(interval);
			}
			else
			{
				intervals.push_back(interval);
			}
		}

		if (!lefts.empty())
		{
			left = std::unique_ptr<intervalTree>(new intervalTree(lefts, depth, minbucket, leftp, centerp));
		}
		if (!rights.empty())
		{
			right = std::unique_ptr<intervalTree>(new intervalTree(rights, depth, minbucket, centerp, rightp));
		}
	}
	*/
}

template <class T, typename K>
IntervalTree<T,K>& 
IntervalTree<T,K>::operator=
(const typename IntervalTree<T,K>::intervalTree& other)
{
	center = other.center;
	intervals = other.intervals;
	left = other.left ? copyTree(*other.left) : nullptr;
	right = other.right ? copyTree(*other.right) : nullptr;
	return *this;
}

template <class T, typename K>
typename IntervalTree<T,K>::intervalVector
IntervalTree<T,K>::findOverlapping
(K start, K stop) const
{
	intervalVector ov;
	this->findOverlapping(start, stop, ov);
	return ov;
}

template <class T, typename K>
void
IntervalTree<T,K>::findOverlapping
(K start, K stop, intervalVector& overlapping) const
{
	if(!intervals.empty() && ! (stop < intervals.front().start))
	{
		for(typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i)
		{
			const interval& interval = *i;
			if (interval.stop >= start && interval.start <= stop)
				overlapping.push_back(interval);
		}
	}

	if(left && start < center)
		left->findOverlapping(start, stop, overlapping);

	if(right && stop > center)
		right->findOverlapping(start, stop, overlapping);

}

template <class T, typename K>
typename IntervalTree<T,K>::intervalVector
IntervalTree<T,K>::findContained
(K start, K stop) const
{
	intervalVector contained;
	this->findContained(start, stop, contained);
	return contained;
}

template <class T, typename K>
void
IntervalTree<T,K>::findContained
(K start, K stop, intervalVector& contained) const
{
	if(!intervals.empty() && ! (stop < intervals.front().start))
	{
		for(typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i)
		{
			const interval& interval = *i;
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
std::unique_ptr<typename IntervalTree<T,K>::intervalTree>
IntervalTree<T,K>::copyTree
(const typename IntervalTree<T,K>::intervalTree& orig)
{
	return std::unique_ptr<intervalTree>(new intervalTree(orig));
}
