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

# ./mobius_exp_mag mob-exp-mag-1K 400 400 5000 0 0 1000
# ../../generate/renorm mob-exp-mag-1K rg 6.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mob-exp-mag-1K.png

../../generate/renorm mobius-12K rg 9.0
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mob-exp-mag-12K.png

# ------------------
# ./divisor_exp_mag divisor-exp-mag-120 400 400 5000 0 0 120
# ../../generate/renorm divisor-exp-mag-120 rg 6.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > div-exp-mag-120.png

# ./divisor_exp_mag divisor-exp-mag-1K 400 400 5000 0 0 1000
# ../../generate/renorm divisor-exp-mag-1K rg 4.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > div-exp-mag-1K.png

# ./divisor_big divisor-12K 400 400 5000 0 0 12000
# ../../generate/renorm divisor-12K rg 1.6
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > div-exp-mag-12K.png

# ../../generate/renorm divisor-120K rg 1.6
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > div-exp-mag-120K.png

# ------------------
# ./liouv_omega_exp_mag liouv-omega-120 400 400 5000 0 0 120
# ../../generate/renorm liouv-omega-120 rg 4.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > liouv-omega-120.png

# ./liouv_omega_exp_mag liouv-omega-1K 400 400 5000 0 0 1000
# ../../generate/renorm liouv-omega-1K rg 4.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > liouv-omega-1K.png

# ------------------
# ./liouv_lambda liouv-lambda-120 400 400 5000 0 0 120
# ../../generate/renorm liouv-lambda-120 rg 4.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > liouv-lambda-120.png

# ./liouv_lambda liouv-lambda-1K 400 400 5000 0 0 1000
# ../../generate/renorm liouv-lambda-1K rg 4.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > liouv-lambda-1K.png

# ------------------
# ./mertens_m mertens-m-120 400 400 5000 0 0 120
# ../../generate/renorm mertens-m-120 rg 4.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mertens-m-120.png

# ./mertens_m mertens-m-1K 400 400 5000 0 0 1000
# ../../generate/renorm mertens-m-1K rg 4.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mertens-m-1K.png

# ------------------
# ./mangoldt_lambda mangoldt-lambda-120 400 400 5000 0 0 120
../../generate/renorm mangoldt-lambda-120 rg 3.0
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mangoldt-lambda-120.png

# ./mangoldt_lambda mangoldt-lambda-1K 400 400 5000 0 0 1000
../../generate/renorm mangoldt-lambda-1K rg 3.0
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mangoldt-lambda-1K.png

# ------------------
# ./thue_morse mertens-m-120 400 400 5000 0 0 120
# ../../generate/renorm thue-morse-120 rg 4.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > thue-morse-120.png

# ./thue_morse thue-morse-1K 400 400 5000 0 0 1000
# ../../generate/renorm thue-morse-1K rg 12.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > thue-morse-1K.png

# ./thue_morse_big thue-morse-1K 400 400 5000 0 0 12000
# ../../generate/renorm thue-morse-12K rg 12.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > thue-morse-12K.png

# ../../generate/renorm thue-morse-120K rg 24.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > thue-morse-120K.png

# ------------------
../../generate/renorm exp-x rg 3.0
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp-x.png

