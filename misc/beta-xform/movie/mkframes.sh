#! /bin/bash

i=1
while [ $i -lt 500 ];
do
	let j=$i
	n=$(echo "0.5 + 0.001 * $j" |bc)
	time ./betadisk beta-egf-$j 800 800 $j 0 0 900 $n &
	let j=$i+1
	n=$(echo "0.5 + 0.001 * $j" |bc)
	time ./betadisk beta-egf-$j 800 800 $j 0 0 900 $n &
	let j=$i+2
	n=$(echo "0.5 + 0.001 * $j" |bc)
	time ./betadisk beta-egf-$j 800 800 $j 0 0 900 $n &
	let j=$i+3
	n=$(echo "0.5 + 0.001 * $j" |bc)
	time ./betadisk beta-egf-$j 800 800 $j 0 0 900 $n &
	let j=$i+4
	n=$(echo "0.5 + 0.001 * $j" |bc)
	time ./betadisk beta-egf-$j 800 800 $j 0 0 900 $n &
	let j=$i+5
	n=$(echo "0.5 + 0.001 * $j" |bc)
	time ./betadisk beta-egf-$j 800 800 $j 0 0 900 $n &
	let j=$i+6
	n=$(echo "0.5 + 0.001 * $j" |bc)
	time ./betadisk beta-egf-$j 800 800 $j 0 0 900 $n &
	let j=$i+7
	n=$(echo "0.5 + 0.001 * $j" |bc)
	time ./betadisk beta-egf-$j 800 800 $j 0 0 900 $n &
	let j=$i+8
	n=$(echo "0.5 + 0.001 * $j" |bc)
	time ./betadisk beta-egf-$j 800 800 $j 0 0 900 $n &
	let i=$i+9
	# sleep 4
done
