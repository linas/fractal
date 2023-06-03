#! /usr/bin/env python3
#
# Solve the equations for a radially infalling flashlight,
# followed by a radially infalling astronaut. Record what
# the coordinates of each. All values are relative to the
# proper time of the flashlight, when the flashlight flashed.
#

import math

# Schwarschild scaled mass
mass = 1

# Initial location of lab frame
rnaught = 20

# Initial location of flashlight
deltar = -0.1

# Proper time step size
taustep = 0.5

print("#\n# Radial infall simulation")
print("# Mass=", mass)
print("# r_0=", rnaught)
print("# delta r_0=", deltar)

# Initial flashlight radial coordinate
rfnaught = rnaught + deltar

# Timelike object radial proper time scale
rscale = 1.5 * math.sqrt(2.0*mass)

# Return the radial coord in the Schwarzschild chart of a radially
# infalling massive object, initially located at `rzero`, at proper
# time `prop`.
def sch_rp_coord(rzero, prop):
	robj = math.pow(rzero, 1.5) - rscale * prop
	if robj < 0:
		return robj
	robj = math.pow(robj, 2.0/3.0)
	return robj

# Return the "time" coord in the Schwarzschild chart of a radially
# infalling massive object, currrently located at the radial
# coordinate `rcoord`.
def sch_tp_coord(rcoord):
	squ = math.sqrt(0.5*rcoord/mass)
	la = (squ - 1.0) / (squ + 1)
	if la < 0.0:
		la = -la
	lg = 2.0*mass* math.log(la)

	time = 2.0 * (rcoord + 6.0*mass) * squ / 3.0 + lg
	return time


tau = 0.0
while True:

	# flashlight radial coordinate
	rflash = sch_rp_coord(rfnaught, tau)
	if rflash < 0.0:
		break

	# flashlight time coord
	tflash = sch_tp_coord(rfnaught) - sch_tp_coord(rflash)

	print("duuude", \
		"{:10.4f}".format(tau), \
		"{:10.4f}".format(rflash), \
		"{:10.4f}".format(tflash) \
		)
	tau += taustep
