#! /bin/bash

#
# ------------------
../../generate/renorm exp-x rg 2.0
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp-x.png

