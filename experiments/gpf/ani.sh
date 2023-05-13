#! /bin/bash

# for i in `seq 0 9` ;
# do
# 	cat gpf-exp-deriv-$i.flo | ../../image/flo2mtv |mtvtoppm | ppmtojpeg > gpf-exp-derive-00$i.jpeg
# done

for i in `seq 10 99` ;
do
	cat gpf-exp-deriv-$i.flo | ../../image/flo2mtv |mtvtoppm | ppmtojpeg > gpf-exp-derive-0$i.jpeg
done

for i in `seq 100 107` ;
do
	cat gpf-exp-deriv-$i.flo | ../../image/flo2mtv |mtvtoppm | ppmtojpeg > gpf-exp-derive-$i.jpeg
done

convert -delay 20 -loop 0 gpf-exp-derive-???.jpeg gpf-exp-deriv-ani.gif
