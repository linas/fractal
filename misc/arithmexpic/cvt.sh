#! /bin/bash

# ./totient_ord_phase tot-ord-phase 400 400 5000 0 0 2
# cat tot-ord-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > tot-ord-phase.png

# ./totient_exp_phase tot-exp-phase 400 400 5000 0 0 60
# cat tot-exp-phase.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > tot-exp-phase.png

# ./totient_exp_mag tot-exp-mag-120 400 400 5000 0 0 120
#../../generate/renorm tot-exp-mag-120 rg 1
#cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > tot-exp-mag-120.png

# ./totient_exp_mag tot-exp-mag-1K 400 400 5000 0 0 1000
# ../../generate/renorm tot-exp-mag-1K rg 0.5
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > tot-exp-mag-1K.png

# ------------------
# ./mobius_exp_mag mob-exp-mag-120 400 400 5000 0 0 120
# ../../generate/renorm mob-exp-mag-120 rg 6.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mob-exp-mag-120.png

# ./mobius_exp_mag mob-exp-mag-120 400 400 5000 0 0 1000
# ../../generate/renorm mob-exp-mag-1K rg 6.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mob-exp-mag-1K.png

# ------------------
# ./divisor_exp_mag div-exp-mag-120 400 400 5000 0 0 120
../../generate/renorm divisor-exp-mag-120 rg 6.0
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > div-exp-mag-120.png

# ./divisor_exp_mag div-exp-mag-120 400 400 5000 0 0 1000
../../generate/renorm divisor-exp-mag-1K rg 4.0
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > div-exp-mag-1K.png

../../generate/renorm exp-x rg 6.0
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp-x.png

