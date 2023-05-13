#! /bin/bash

cat revsum.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > revsum.png
