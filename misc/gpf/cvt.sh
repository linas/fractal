#! /bin/bash

# cat gpf-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-phase.png

# ../../generate/takelog gpf-real log-real
# ../../generate/renorm log-real rg 0.24 1
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-real.png

# ../../generate/takelog gpf-abs log-abs
# cat log-abs.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-abs.png

# cat gpf-exp-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-phase.png
# cat gpf-exp-phase-fixed.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-phase-fixed.png

# ../../generate/takelog gpf-exp-abs log-exp-abs
# ../../generate/renorm log-exp-abs rg 0.0124
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs.png

# ../../generate/renorm gpf-exp-abs rg 0.024
# ../../generate/renorm gpf-exp-abs rg 0.24
# ../../generate/renorm gpf-exp-abs rg 0.4
# ../../generate/renorm gpf-exp-abs rg 0.9999999
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs.png
# ../../generate/renorm gpf-exp-abs-fixed rg 0.24
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs-0.24-fixed.png

# ../../generate/renorm gpf-exp-x rg 480
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-x.png

# ../../generate/renorm gpf-exp-abs-big-fixed rg 0.4
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-ex-abs-big-fixed.png

# huge == width of 4320
# ../../generate/renorm gpf-exp-abs-huge rg 0.03124
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs-huge.png

# ../../generate/renorm gpf-exp-abs-4320-fixed rg 0.03124
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs-4320-fixed.png

# ../../generate/renorm gpf-exp-abs-40K rg 0.0041124
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs-40K.png

# ../../generate/renorm gpf-exp-abs-40K-fixed rg 0.0041124
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs-40K-fixed.png

# ../../generate/renorm gpf-exp-abs-250K-uni rg 0.4
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs-250K-uni.png

# ../../generate/renorm rand-exp-fixed rg 1.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > rand-exp-fixed.png

# ../../generate/renorm rand-exp-prime-fixed rg 1.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > rand-exp-prime-fixed.png
# ../../generate/renorm gpf-exp-x rg 4.0
# ../../generate/renorm gpf-exp-x rg 1.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-x.png

# ../../generate/renorm gpf-exp-x rg 0.000000004
# ../../generate/renorm recip-flat-13 rg 24.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-flat-13.png
# ../../generate/renorm recip-flat-16 rg 24.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-flat-16.png

# ../../generate/renorm recip-flat-16-803 rg 24.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-flat-16-803.png

# ../../generate/renorm recip-flat-10-16 rg 24.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-flat-10-16.png

# ../../generate/renorm recip-flat-wide-15 rg 24.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-flat-wide-15.png

# enscribe some lines ---------
# ../../generate/renorm recip-flat-16 rg 24.0
# cat hyper.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > hyper.png
# ../../generate/mask rg hyper rma
# cat rma.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-flat-16-scribed.png

# cat recip.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip.png

#../../generate/renorm recip-exp rg 0.24
#cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-exp.png


# ../../generate/renorm recip-exp-4320 rg 0.5
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-exp-4320.png
../../generate/renorm recip-exp-4320-fixed rg 0.3333
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-exp-4320-fixed.png

# ../../generate/renorm recip-exp-40K rg 3.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-exp-40K.png

# ../../generate/renorm recip-exp-40K-fixed rg 0.5
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-exp-40K-fixed.png

# ../../generate/renorm recip-exp-41K rg 24.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-exp-41K.png

# ../../generate/renorm recip-exp-240K rg 16.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-exp-240K.png

# ../../generate/renorm exp-close rg 0.5
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp-close.png

# ../../generate/renorm diri rg 100.0 -0.5
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > diri-x.png

# ../../generate/renorm poch-rise rg 0.2
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > poch-rise.png

# ../../generate/renorm poch-rise-321M rg 0.2
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > poch-rise-321M.png

# ../../generate/renorm poch-rise-x rg 1
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > poch-rise-x.png

# ../../generate/renorm poch-fall rg 0.2
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > poch-fall.png

# ../../generate/renorm exp+0.5 rg 0.000025
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp+0.5.png

# ../../generate/renorm exp+1.0 rg 0.44
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp+1.0.png

# ../../generate/renorm exp-re rg 0.4
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp-re.png

# ../../generate/renorm exp-re-2K rg 0.4
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp-re-2K.png

# ../../generate/renorm exp-re-2160 rg 0.4
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp-re-2160.png

# ../../generate/renorm exp-re-4K rg 0.4
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp-re-4K.png

# ../../generate/renorm exp-re-8K rg 0.4
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp-re-8K.png

# ../../generate/renorm exp-2.0 rg 0.3
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp-2.0.png

# ../../generate/renorm x rg 1
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > x.png

# ../../generate/renorm gpf-exp-abs-x rg 1
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs-x.png

# cat gpf-norm-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-norm-phase.png
# cat gpf-pois-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-pois-phase.png
# cat gpf-lamb-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-lamb-phase.png
