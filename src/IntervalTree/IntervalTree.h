#ifndef __INTERVAL_TREE_H
#define __INTERVAL_TREE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>

#include <boost/serialization/serialization.hpp>

using namespace std;
using namespace boost;


template <class T>
class Interval
{
public:
	int start;
	int stop;
	T value;
	Interval(int s, int e, const T& v) : start(s), stop(e), value(v) {}
	Interval() : start(0), stop(0) {}
	
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & start;
		ar & stop;
		ar & value;
	}
};

template <class T>
bool StartCmp(const Interval<T>& a, const Interval<T>& b)
{
	return a.start < b.start;
}

template <class T>
bool StopCmp(const Interval<T>& a, const Interval<T>& b)
{
	return a.stop < b.stop;
}

template <class T>
ostream& operator<<(ostream& out, const Interval<T>& i)
{
	out << "Interval(" << i.start << ", " << i.stop << "): " << i.value;
	return out;
}

template <class T>
ostream& operator<<(ostream& out, const vector<Interval<T> >& intervals)
{
	for (typename vector<Interval<T> >::const_iterator intervalIter = intervals.begin(); intervalIter != intervals.end(); intervalIter++)
	{
		out << *intervalIter << endl;
	}
	return out;
}

template <class T>
class NearestAccumulator
{
public:
	NearestAccumulator() : mNearestDistance(numeric_limits<int>::max()) {}
										
	void Add(int distance, const T& value)
	{
		if (distance < mNearestDistance)
		{
			mNearestDistance = distance;
			mNearestValues.clear();
		}
		
		if (distance <= mNearestDistance)
		{
			mNearestValues.push_back(value);
		}
	}
	
	const vector<T>& Get()
	{
		return mNearestValues;
	}
	
private:
	int mNearestDistance;
	vector<T> mNearestValues;
};

template <class T>
class IntervalTree
{
public:
	IntervalTree<T>() : mLeft(0), mRight(0), mCenter(0) {}
	
	IntervalTree<T>(const IntervalTree<T>& other)
	{
		mCenter = other.mCenter;
		
		mIntervals = other.mIntervals;
		
		if (other.mLeft)
		{
			mLeft = new IntervalTree<T>(*other.mLeft);
		}
		else
		{
			mLeft = 0;
		}
		
		if (other.mRight)
		{
			mRight = new IntervalTree<T>(*other.mRight);
		}
		else
		{
			mRight = 0;
		}
	}
	
	IntervalTree<T>& operator=(const IntervalTree<T>& other)
	{
		mCenter = other.mCenter;
		
		mIntervals = other.mIntervals;
		
		delete mLeft;
		if (other.mLeft)
		{
			mLeft = new IntervalTree<T>(*other.mLeft);
		}
		else
		{
			mLeft = 0;
		}
		
		delete mRight;
		if (other.mRight)
		{
			mRight = new IntervalTree<T>(*other.mRight);
		}
		else
		{
			mRight = 0;
		}
		
		return *this;
	}
	
	IntervalTree<T>(vector<Interval<T> >& intervals, 
					unsigned int maxdepth = 16,
					unsigned int minbucket = 64,
					unsigned int maxbucket = 512)
	: mLeft(0), mRight(0)
	{
		sort(intervals.begin(), intervals.end(), StartCmp<T>);
		
		int leftextent = intervals.front().start;
		int rightextent = max_element(intervals.begin(), intervals.end(), StopCmp<T>)->stop;
		
		Construct(intervals, maxdepth, minbucket, maxbucket, leftextent, rightextent);
	}
	
	void FindOverlapping(int start, int stop, vector<T>& overlapping) const
	{
		if (stop >= mIntervals.front().start)
		{
			for (typename vector<Interval<T> >::const_iterator intervalIter = mIntervals.begin(); intervalIter != mIntervals.end(); intervalIter++)
			{
				if (intervalIter->stop >= start && intervalIter->start <= stop)
				{
					overlapping.push_back(intervalIter->value);
				}
			}
		}
		
		if (mLeft && start <= mCenter)
		{
			mLeft->FindOverlapping(start, stop, overlapping);
		}
		
		if (mRight && stop >= mCenter)
		{
			mRight->FindOverlapping(start, stop, overlapping);
		}
	}
	
	void FindContained(int start, int stop, vector<T>& contained) const
	{
		if (stop >= mIntervals.front().start)
		{
			for (typename vector<Interval<T> >::const_iterator intervalIter = mIntervals.begin(); intervalIter != mIntervals.end(); intervalIter++)
			{
				if (intervalIter->start >= start && intervalIter->stop <= stop)
				{
					contained.push_back(intervalIter->value);
				}
			}
		}
		
		if (mLeft && start <= mCenter)
		{
			mLeft->FindContained(start, stop, contained);
		}
		
		if (mRight && stop >= mCenter)
		{
			mRight->FindContained(start, stop, contained);
		}
		
	}
	
	void FindNearest(int position, vector<T>& nearest) const
	{
		NearestAccumulator<T> nearestAccumulator;
		FindNearest(position, nearestAccumulator);
		nearest.insert(nearest.end(), nearestAccumulator.Get().begin(), nearestAccumulator.Get().end());
	}
	
	~IntervalTree()
	{
		delete mLeft;
		delete mRight;
	}
	
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & mIntervals;
		ar & mLeft;
		ar & mRight;
		ar & mCenter;
	}
	
private:
	IntervalTree<T>(vector<Interval<T> >& intervals, 
					unsigned int maxdepth,
					unsigned int minbucket,
					unsigned int maxbucket,
					int leftextent,
					int rightextent)
	: mLeft(0), mRight(0)
	{
		Construct(intervals, maxdepth, minbucket, maxbucket, leftextent, rightextent);
	}
	
	void Construct(vector<Interval<T> >& intervals, unsigned int depth, unsigned int minbucket,
				   unsigned int maxbucket, int leftextent, int rightextent)
	{
		if (depth == 1 || (intervals.size() < minbucket && intervals.size() < maxbucket))
		{
			mIntervals = intervals;
		}
		else
		{
			mCenter = intervals.at(intervals.size() / 2).start;
			
			vector<Interval<T> > lefts;
			vector<Interval<T> > rights;
			
			for (typename vector<Interval<T> >::iterator i = intervals.begin(); i != intervals.end(); ++i)
			{
				Interval<T>& interval = *i;
				if (interval.stop < mCenter)
				{
					lefts.push_back(interval);
				}
				else if (interval.start > mCenter)
				{
					rights.push_back(interval);
				}
				else
				{
					mIntervals.push_back(interval);
				}
			}
			
			if (!lefts.empty())
			{
				mLeft = new IntervalTree<T>(lefts, depth - 1, minbucket, maxbucket, leftextent, mCenter);
			}
			if (!rights.empty())
			{
				mRight = new IntervalTree<T>(rights, depth - 1, minbucket, maxbucket, mCenter, rightextent);
			}
		}
	}
	
	void FindNearest(int position, NearestAccumulator<T>& nearest) const
	{
		for (typename vector<Interval<T> >::const_iterator intervalIter = mIntervals.begin(); intervalIter != mIntervals.end(); intervalIter++)
		{
			if (position < intervalIter->start)
			{
				nearest.Add(intervalIter->start - position, intervalIter->value);
			}
			else if (position > intervalIter->stop)
			{
				nearest.Add(position - intervalIter->stop, intervalIter->value);
			}
			else
			{
				nearest.Add(0, intervalIter->value);
			}
		}
		
		if (mLeft && position < mCenter)
		{
			mLeft->FindNearest(position, nearest);
		}
		
		if (mRight && position > mCenter)
		{
			mRight->FindNearest(position, nearest);
		}
	}
	
	vector<Interval<T> > mIntervals;
	IntervalTree<T>* mLeft;
	IntervalTree<T>* mRight;
	int mCenter;
};

#endif

