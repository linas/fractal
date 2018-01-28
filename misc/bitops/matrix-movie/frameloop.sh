#! /bin/bash


for i in `seq 1 3`;
do
	# let K=0.6+0.1*$i
	K="$(bc <<<  "scale=3; 0.6 + 0.1 * $i")"
	echo "do frame $i K= $K"
	../matrix foo 200 200 1 0 100.0 200 $K
	../../../generate/renorm foo j 1.5 0
	cat j.flo | ../../../image/flo2mtv |mtvtoppm | pnmtopng > matrix-$K.png
done    


