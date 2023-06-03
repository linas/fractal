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
rfncube = math.pow(rfnaught, 1.5)

# Timelike object radial proper time scale
rscale = 1.5 * math.sqrt(2.0*mass)

# Return the radial coord in the Schwarzschild chart of an infalling
# massive object, initially located at `rzero`, at proper time `prop`.
def sch_rp_coord(rzero, prop):
	robj = math.pow(rzero, 1.5) - rscale * prop
	if robj < 0:
		return robj
	robj = math.pow(robj, 2.0/3.0)
	return robj


tau = 0.0
while True:

	rflash = rfncube - rscale * tau
	if rflash < 0:
		break;
	rflash = math.pow(rflash, 2.0/3.0)

	print("duuude", \
		"{:10.4f}".format(tau), \
		"{:10.4f}".format(rflash) \
		)
	tau += taustep
