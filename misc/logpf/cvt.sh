#! /bin/bash

#
# ------------------
../../generate/renorm exp-x rg 0.025
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp-x.png

