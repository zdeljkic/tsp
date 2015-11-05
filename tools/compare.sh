#!/bin/bash

TESTGEN_BIN=../testgen/bin/Debug/testgen
TSP_BIN=../tsp/bin/Debug/tsp

testTSP()
{
	total_dist=0
	total_time=0
	for i in $(seq 1 $1)
	do
		dist=$({ time $TSP_BIN $2 <testIn$i.txt; } 2>time.txt | tail -n 1 | cut -d ' ' -f 2)
		time_str=$(cat time.txt | tail -n 3 | head -n 1 | cut -f 2)
		min=$(echo $time_str | cut -d m -f 1)
		sec=$(echo $time_str | cut -d m -f 2 | cut -d s -f 1)
		time=$(echo "60*$min + $sec" | bc)
		total_time=$(echo "$total_time + $time" | bc)
		total_dist=$(echo "$total_dist + $dist" | bc)
	done
	
	avg_dist=$(echo "scale=3; $total_dist / $1" | bc)
	avg_time=$(echo "scale=3; $total_time / $1" | bc)
	
	echo "Average distance: $avg_dist"
	echo "Average time: ${avg_time}s"
	
	rm time.txt
}

if [ $# != 2 ]
then
	echo 'Usage: $0 number_of_test_samples number_of_test_repetitions'
	exit 1
fi

echo "Generating $2 test inputs of size $1..."
for i in $(seq 1 $2)
do
	$TESTGEN_BIN $1 > testIn$i.txt
done

echo
echo 'Testing nearest neighbor:'
testTSP $2 n

echo
echo 'Testing greedy:'
testTSP $2 g

echo
echo 'Testing greedy + 2opt:'
testTSP $2 2

rm testIn*.txt
