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

# cat gpf-norm-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-norm-phase.png
# cat gpf-pois-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-pois-phase.png
cat gpf-lamb-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-lamb-phase.png
