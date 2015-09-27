#! /bin/bash

# cat lat.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > lat.png
# Below is OK for 10
# ../../generate/renorm lat rlat 0.08
# Below is nice for 12
# ../../generate/renorm lat rlat 0.02
# ../../generate/renorm lat rlat 0.024
# ../../generate/renorm olat rlat 0.024
../../generate/renorm olat rlat 0.001
# Below is nice for 16
# ../../generate/renorm lat rlat 0.002
cat rlat.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > olat.png


# ../../generate/recip olat clat
# ../../generate/takelog clat llat
# ../../generate/renorm llat rlat -0.15
# cat rlat.flo | ../../image/flo2mtv |mtvtoppm | pnmtopng > olat.png
