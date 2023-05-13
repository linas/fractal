#! /bin/bash

for i in `seq 1 9` ;
do
	cat beta-egf-$i.flo | ../../../image/flo2mtv |mtvtoppm | ppmtojpeg > beta-egf-00$i.jpeg
done

for i in `seq 10 99` ;
do
	cat beta-egf-$i.flo | ../../../image/flo2mtv |mtvtoppm | ppmtojpeg > beta-egf-0$i.jpeg
done

for i in `seq 100 504` ;
do
	cat beta-egf-$i.flo | ../../../image/flo2mtv |mtvtoppm | ppmtojpeg > beta-egf-$i.jpeg
done

convert -delay 20 -loop 0 beta-egf-???.jpeg beta-egf-ani.gif
