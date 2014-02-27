import unittest
import pygenes
import random
import itertools

def overlapping(interval, start, stop):
	if interval[2] < start or interval[1] > stop:
		return False
	return True

def contained(interval, start, stop):
	if interval[1] >= start and interval[2] <= stop:
		return True
	return False
	
def distance(interval, position):
	if position < interval[1]:
		return interval[1] - position
	elif position > interval[2]:
		return position - interval[2]
	else:
		return 0

def find_overlapping(intervals, start, stop):
	for interval in intervals:
		if overlapping(interval, start, stop):
			yield interval[0]

def find_contained(intervals, start, stop):
	for interval in intervals:
		if contained(interval, start, stop):
			yield interval[0]

def find_nearest(intervals, position):
	sorted_intervals = sorted(intervals, key=lambda interval: distance(interval, position))
	nearest_intervals = next(itertools.groupby(sorted_intervals, key=lambda interval: distance(interval, position)))[1]
	for interval in nearest_intervals:
		yield interval[0]

class intervaltree_test(unittest.TestCase):
	
	def setUp(self):
		pass
	
	def test_simple(self):
		
		intervals = list()
		intervals.append((1, 3, 5))
		intervals.append((2, 9, 10))
		intervals.append((3, 10, 14))
		interval_tree = pygenes.IntervalTree(intervals)
		
		overlapping = interval_tree.find_overlapping(6, 12)
		contained = interval_tree.find_contained(6, 12)
		nearest = interval_tree.find_nearest(7)
		
		overlapping_ids = set([id for id in overlapping])
		contained_ids = set([id for id in contained])
		nearest_ids = set([id for id in nearest])
		
		self.assertEqual(overlapping_ids, set([2, 3]))
		self.assertEqual(contained_ids, set([2]))
		self.assertEqual(nearest_ids, set([1, 2]))
		
	def construct_random_intervals(self, num_intervals, max_start, max_size):
		intervals = list()
		for idx in range(num_intervals):
			start = random.randint(0, max_start)
			stop = start + random.randint(0, max_size - 1)
			intervals.append((idx, start, stop))
		return intervals
		
	def compare_random_intervals(self, num_intervals, max_start, max_size, intervals, interval_tree_func, test_func):
		for idx in range(num_intervals):
			start = random.randint(0, max_start)
			stop = start + random.randint(0, max_size - 1)
			query_ids = set([id for id in interval_tree_func(start, stop)])
			self.assertEqual(query_ids, set(test_func(intervals, start, stop)))
	
	def compare_random_positions(self, num_positions, max_position, intervals, interval_tree_func, test_func):
		for idx in range(num_positions):
			position = random.randint(0, max_position)
			query_ids = set([id for id in interval_tree_func(position)])
			self.assertEqual(query_ids, set(test_func(intervals, position)))
	
	def test_random_dense(self):
		intervals = self.construct_random_intervals(1000, 1000, 100)
		interval_tree = pygenes.IntervalTree(intervals)
		self.compare_random_intervals(100, 10000, 100, intervals, interval_tree.find_overlapping, find_overlapping)
		self.compare_random_intervals(100, 10000, 100, intervals, interval_tree.find_contained, find_contained)
		self.compare_random_positions(100, 10000, intervals, interval_tree.find_nearest, find_nearest)
		
	def test_random_sparse(self):
		intervals = self.construct_random_intervals(1000, 10000, 10)
		interval_tree = pygenes.IntervalTree(intervals)
		self.compare_random_intervals(100, 10000, 100, intervals, interval_tree.find_overlapping, find_overlapping)
		self.compare_random_intervals(100, 10000, 100, intervals, interval_tree.find_contained, find_contained)
		self.compare_random_positions(100, 10000, intervals, interval_tree.find_nearest, find_nearest)
		
		
if __name__ == '__main__':
	unittest.main()
	
