#! /bin/bash

# cat gpf-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-phase.png

# ../../generate/takelog gpf-real log-real
# ../../generate/renorm log-real rg 0.24 1
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-real.png

../../generate/takelog gpf-abs log-abs
cat log-abs.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-abs.png
