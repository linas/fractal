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

# ../../generate/renorm totient-12K rg 0.0667
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > tot-exp-mag-12K.png

# ./totient_big totient-log-12K 200 200 5000 0 0 12000
# ../../generate/renorm totient-log-12K rg 0.5
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > totient-log-12K.png

# ../../generate/renorm totient-120K rg 0.125
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > totient-120K.png

# ------------------
# ./carmichael carmichael-120 400 400 50005 0 0 120
# ../../generate/renorm carmichael-120 rg 0.25
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > carmichael-120.png

# ../../generate/renorm carmichael-1K rg 0.12
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > carmichael-1K.png

# ../../generate/renorm carmichael-12K rg 0.5
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > carmichael-12K.png

# ../../generate/renorm carmichael-120K rg 0.18
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > carmichael-120K.png

# ------------------
# ./mobius_exp_mag mob-exp-mag-120 400 400 5000 0 0 120
# ../../generate/renorm mob-exp-mag-120 rg 6.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mob-exp-mag-120.png

# ./mobius_exp_mag mob-exp-mag-1K 400 400 5000 0 0 1000
# ../../generate/renorm mob-exp-mag-1K rg 6.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mob-exp-mag-1K.png

# ../../generate/renorm mobius-12K rg 9.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mob-exp-mag-12K.png

# ../../generate/renorm mobius-120K rg 9.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mob-exp-mag-120K.png

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

# ../../generate/renorm sigma-one-120 rg 2.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > sigma-one-120.png

# ../../generate/renorm sigma-one-1K rg 0.8
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > sigma-one-1K.png

# ../../generate/renorm sigma-two-120 rg 1.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > sigma-two-120.png

# ../../generate/renorm sigma-two-1K rg 0.16
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > sigma-two-1K.png

# ------------------
# ./little_omega little-omega-120 400 400 5000 0 0 120
# ../../generate/renorm little-omega-120 rg 4.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > little-omega-120.png
#
# ../../generate/renorm little-omega-1K rg 8.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > little-omega-1K.png
#
# ../../generate/renorm little-omega-12K rg 16.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > little-omega-12K.png
#
# ../../generate/renorm little-omega-120K rg 16.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > little-omega-120K.png

# ------------------
# ./unitary unitary-120 400 400 5000 0 0 120
# ../../generate/renorm unitary-120 rg 2.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > unitary-120.png
#
# ../../generate/renorm unitary-1K rg 2.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > unitary-1K.png
#
# ../../generate/renorm unitary-12K rg 3.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > unitary-12K.png
#
# ../../generate/renorm unitary-120K rg 3.0
#  cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > unitary-120K.png

# ------------------
# ./big_omega_exp_mag liouv-omega-120 400 400 5000 0 0 120
# ../../generate/renorm liouv-omega-120 rg 4.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > liouv-omega-120.png

# ./big_omega_exp_mag liouv-omega-1K 400 400 5000 0 0 1000
# ../../generate/renorm liouv-omega-1K rg 4.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > liouv-omega-1K.png

# ../../generate/renorm liouv-omega-12K rg 8.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > liouv-omega-12K.png

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
# ../../generate/renorm mangoldt-lambda-120 rg 3.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mangoldt-lambda-120.png

# ./mangoldt_lambda mangoldt-lambda-1K 400 400 5000 0 0 1000
# ../../generate/renorm mangoldt-lambda-1K rg 3.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mangoldt-lambda-1K.png

# ../../generate/renorm mangoldt-lambda-12K rg 3.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > mangoldt-lambda-12K.png

# ------------------
# ./exp_mangoldt_lambda emango-lambda-120 400 400 5000 0 0 120
# ../../generate/renorm emango-lambda-120 rg 2
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > emango-lambda-120.png

# ./exp_emango_lambda emango-lambda-1K 400 400 5000 0 0 1000
# ../../generate/renorm emango-lambda-1K rg 1.2
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > emango-lambda-1K.png

# ../../generate/renorm emango-lambda-12K rg 0.4
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > emango-lambda-12K.png

# ../../generate/renorm emango-lambda-120K rg 0.1
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > emango-lambda-120K.png

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
# ./isqrt_bit isqrt-120 400 400 5000 0 0 120
# ../../generate/renorm isqrt-120 rg 8.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > isqrt-120.png
#
# ../../generate/renorm isqrt-1K rg 16.0
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > isqrt-1K.png
#
# ------------------
# ./partition_bit partition-120 400 400 5000 0 0 120
# ../../generate/renorm partition-12 rg 1
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > partition-12.png

# ../../generate/takeroot partition-120 rg 6
# ../../generate/renorm rg rg 0.8
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > partition-120.png

# ../../generate/renorm partition-160 rg 1
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > partition-160.png

# ../../generate/renorm partition-1200 rg 1
# cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > partition-1200.png

../../generate/takeroot partition-1K rg 6
../../generate/renorm rg rg 0.8
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > partition-1K.png
#
# ------------------
../../generate/renorm exp-x rg 3.0
cat rg.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > exp-x.png

