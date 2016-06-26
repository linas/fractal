#! /bin/bash

i=0
while [ $i -lt 100 ];
do
	j=$i
	time ./gpf-2d gpf-exp-deriv-$j 600 600 $j 0 0 200 &
	j=$i+1
	time ./gpf-2d gpf-exp-deriv-$j 600 600 $j 0 0 200 &
	j=$i+2
	time ./gpf-2d gpf-exp-deriv-$j 600 600 $j 0 0 200 &
	j=$i+3
	time ./gpf-2d gpf-exp-deriv-$j 600 600 $j 0 0 200 &
	j=$i+4
	time ./gpf-2d gpf-exp-deriv-$j 600 600 $j 0 0 200 &
	j=$i+5
	time ./gpf-2d gpf-exp-deriv-$j 600 600 $j 0 0 200 &
	j=$i+6
	time ./gpf-2d gpf-exp-deriv-$j 600 600 $j 0 0 200 &
	j=$i+7
	time ./gpf-2d gpf-exp-deriv-$j 600 600 $j 0 0 200 &
	j=$i+8
	time ./gpf-2d gpf-exp-deriv-$j 600 600 $j 0 0 200
	i=$i+9
	sleep 6
done
