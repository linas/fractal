#! /bin/bash

cat sreal.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > sreal.png
