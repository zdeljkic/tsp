#!/bin/bash

TESTGEN_BIN=../testgen/bin/Debug/testgen
TSP_BIN=../tsp/bin/Debug/tsp

testTSP()
{
	min_dist=2000000000 # 1000 cities * 2000000 is more than max distance
	max_dist=0
	total_dist=0
	min_time=1000.0
	max_time=0.0
	total_time=0

	for i in $(seq 1 $1)
	do
		dist=$({ time $TSP_BIN $2 <testIn$i.txt; } 2>time.txt | tail -n 1 | cut -d ' ' -f 2)
		time_str=$(cat time.txt | tail -n 3 | head -n 1 | cut -f 2)
		min=$(echo $time_str | cut -d m -f 1)
		sec=$(echo $time_str | cut -d m -f 2 | cut -d s -f 1)
		time=$(echo "60*$min + $sec" | bc)

		min_dist=$(echo "if ($dist < $min_dist) $dist else $min_dist" | bc)
		max_dist=$(echo "if ($dist > $max_dist) $dist else $max_dist" | bc)
		total_dist=$(echo "$total_dist + $dist" | bc)

		min_time=$(echo "scale=3; if ($time < $min_time) $time else $min_time" | bc)
		max_time=$(echo "scale=3; if ($time > $max_time) $time else $max_time" | bc)
		total_time=$(echo "$total_time + $time" | bc)
	done

	avg_dist=$(echo "scale=3; $total_dist / $1" | bc)
	avg_time=$(echo "scale=3; $total_time / $1" | bc)

	echo "min/max/avg distance: $min_dist/$max_dist/$avg_dist"
	echo "min/max/avg time: $min_time/$max_time/$avg_time"

	rm time.txt
}

if [ $# != 2 -a $# != 3 ]
then
	echo 'Usage: $0 number_of_test_samples number_of_test_repetitions [old]'
	exit 1
fi

if [ $# != 3 ]
then
	echo "Generating $2 test inputs of size $1..."
	for i in $(seq 1 $2)
	do
		$TESTGEN_BIN $1 > testIn$i.txt
	done
fi

echo
echo 'Testing nearest neighbor:'
testTSP $2 n

echo
echo 'Testing greedy:'
testTSP $2 g

echo
echo 'Testing randomized greedy:'
testTSP $2 G

echo
echo 'Testing greedy + 2opt:'
testTSP $2 1

echo
echo 'Testing randomized greedy + 2opt:'
testTSP $2 2

echo
echo 'Testing iterative randomized greedy + 2opt:'
testTSP $2 3

if [ $# != 3 ]
then
	rm testIn*.txt
fi
