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
../../generate/renorm gpf-exp-abs-huge rg 0.03124
# ../../generate/renorm gpf-exp-abs-huge rg 0.01124
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-exp-abs-huge.png

# cat gpf-norm-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-norm-phase.png
# cat gpf-pois-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-pois-phase.png
# cat gpf-lamb-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-lamb-phase.png
