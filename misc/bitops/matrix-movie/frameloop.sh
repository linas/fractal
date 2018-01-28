#! /bin/bash


for i in `seq 1 250`;
do
	K="$(bc <<<  "scale=5; 0.5 + $i / 503")"
	echo "do frame $i K= $K"
	../matrix foo 600 600 1 0 300.0 600 $K
	../../../generate/renorm foo j 6.5 0
	../../../generate/abs j j
	cat j.flo | ../../../image/flo2mtv |mtvtoppm | pnmtopng > matrix-0$K.png
done    


