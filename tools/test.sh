#!/bin/bash

TESTGEN_BIN=../testgen/bin/Debug/testgen
TSP_BIN=../tsp/bin/Debug/tsp

testTSP()
{
	total_dist=0
	max_time=0.0

	for i in $(seq 1 $1)
	do
		dist=$({ time $TSP_BIN $2 <${k}testIn${i}.txt; } 2>time.txt | tail -n 1 | cut -d ' ' -f 2)
		time_str=$(cat time.txt | tail -n 3 | head -n 1 | cut -f 2)
		min=$(echo $time_str | cut -d m -f 1)
		sec=$(echo $time_str | cut -d m -f 2 | cut -d s -f 1)
		time=$(echo "60*$min + $sec" | bc)

		total_dist=$(echo "$total_dist + $dist" | bc)
		max_time=$(echo "scale=3; if ($time > $max_time) $time else $max_time" | bc)
	done

	avg_dist=$(echo "scale=3; $total_dist / $1" | bc)

	#echo "$max_time	$avg_dist"

	rm time.txt
}

if [ $# != 2 -a $# != 3 ]
then
	echo 'Usage: $0 step_of_n repeat'
	exit 1
fi

score=0
n=0
for k in $(seq $1 $1 1000)
do
	if [ $# != 3 ]
	then
		for j in $(seq 1 $2)
		do
			$TESTGEN_BIN $k > ${k}testIn${j}.txt
		done
	fi

	testTSP $2 g
	gd=$avg_dist
	testTSP $2 3
	bv=$avg_dist
	
	n=$(($n + 1))
	
	s=$(echo "scale=3; 1 - $bv / $gd" | bc)
	echo "size=$k max_time=$max_time score=$s"
	score=$(echo "scale=3; $score + $s" | bc)
done

score=$(echo "scale=3; $score / $n" | bc)

echo "Average score: $score"

if [ $# != 3 ]
then
	rm *testIn*.txt
fi
