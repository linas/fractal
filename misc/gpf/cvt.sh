#! /bin/bash

# cat gpf-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-phase.png

# ../../generate/takelog gpf-real log-real
# ../../generate/renorm log-real rg 0.24 1
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-real.png

# ../../generate/takelog gpf-abs log-abs
# cat log-abs.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-abs.png

# cat gpf-exp-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-phase.png

# ../../generate/takelog gpf-exp-abs log-exp-abs
# ../../generate/renorm log-exp-abs rg 0.0124
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs.png

# ../../generate/renorm gpf-exp-abs rg 0.024
# ../../generate/renorm gpf-exp-abs rg 0.24
# ../../generate/renorm gpf-exp-abs rg 0.4
# ../../generate/renorm gpf-exp-abs rg 0.9999999
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs.png

# big == width of 720
# ../../generate/renorm gpf-exp-abs-big rg 0.4
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs-big.png
# cat gpf-exp-abs-big.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs-big.png

# huge == width of 4320
# ../../generate/renorm gpf-exp-abs-huge rg 0.03124
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs-huge.png

# ../../generate/renorm gpf-exp-abs-40K rg 0.0041124
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs-40K.png

# ../../generate/renorm gpf-exp-abs-250K-uni rg 0.4
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs-250K-uni.png
# cat gpf-exp-x.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-x.png
# cat recip.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip.png

#../../generate/renorm recip-exp rg 0.24
#cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-exp.png


# ../../generate/renorm recip-exp-4320 rg 0.5
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-exp-4320.png

# ../../generate/renorm recip-exp-40K rg 3.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-exp-40K.png

# ../../generate/renorm recip-exp-41K rg 24.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > recip-exp-41K.png

# ../../generate/renorm exp-close rg 0.5
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp-close.png

# ../../generate/renorm exp+0.5 rg 0.000025
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp+0.5.png

../../generate/renorm exp+1.0 rg 0.4
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp+1.0.png

# ../../generate/renorm x rg 1
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > x.png

# ../../generate/renorm gpf-exp-abs-x rg 1
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs-x.png

# cat gpf-norm-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-norm-phase.png
# cat gpf-pois-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-pois-phase.png
# cat gpf-lamb-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-lamb-phase.png
