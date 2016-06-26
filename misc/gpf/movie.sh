#! /bin/bash

for i in `seq 0 50` ;
do
	time ./gpf-2d gpf-exp-deriv-$i 600 600 $i 0 0 200
done
