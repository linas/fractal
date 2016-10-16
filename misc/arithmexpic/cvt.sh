#! /bin/bash

# cat tot-ord-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > tot-ord-phase.png
# cat tot-exp-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > tot-exp-phase.png

../../generate/renorm tot-exp-mag rg 0.5
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > tot-exp-mag.png

