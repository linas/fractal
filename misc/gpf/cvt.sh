#! /bin/bash

cat gpf-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > gpf-phase.png
