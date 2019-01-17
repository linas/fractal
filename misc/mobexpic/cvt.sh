#! /bin/bash

# Mobius function
# ./dirichlet mob 800 800 5000 0 0 440 -1
# ../../generate/renorm mob rg 4
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mob.png

# Plain-old exponential
# ./dirichlet exp 800 800 5000 0 0 440 1
# ../../generate/renorm exp rg 1
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp.png

# divisor function
# ./dirichlet div 800 800 5000 0 0 440 2
# ../../generate/renorm div rg 2
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > div.png

# unit cubed function
# time ./dirichlet unit3 800 800 5000 0 0 440 3
../../generate/renorm unit3 rg 0.5
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > unit3.png

# Mobius-suqared function
# time ./dirichlet mob2 800 800 5000 0 0 440 -2
# ../../generate/renorm mob2 rg 1
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mob2.png

