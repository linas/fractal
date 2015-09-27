#! /bin/bash

cat lat.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > lat.png
