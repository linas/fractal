#! /bin/bash

#
# ------------------
../../generate/takelog coarse lg
../../generate/renorm lg rg 0.11
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > cor.png

