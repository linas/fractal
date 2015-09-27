#! /bin/bash

# cat lat.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > lat.png
../../generate/renorm lat rlat 0.1
cat rlat.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > lat.png
