#! /usr/bin/env python3
#
# Solve the equations for a radially infalling flashlight,
# followed by a radially infalling astronaut. Record what
# the coordinates of each. All values are relative to the
# proper time of the flashlight, when the flashlight flashed.
#

import math
import scipy
import scipy.optimize

# Schwarschild scaled mass
mass = 1

# Initial location of lab frame
rnaught = 20

# Initial location of flashlight
deltar = -0.1

# Proper time step size
taustep = 0.5

print("#\n# Radial infall simulation")
print("# Mass =", mass)
print("# r_0 =", rnaught)
print("# delta r_0 =", deltar)
print("#")
print("# Columns, in order:")
print("# proper time")
print("# flashlight r")
print("# flashlight t")
print("# astronaut r")
print("# astronaut t")
print("#")

# Initial flashlight radial coordinate
rfnaught = rnaught + deltar

# Timelike object radial proper time scale
rscale = 1.5 * math.sqrt(2.0*mass)

# Return the radial coord in the Schwarzschild chart of a radially
# infalling massive object, initially located at `rzero`, at proper
# time `prop`.
def schw_rp_coord(rzero, prop):
	robj = math.pow(rzero, 1.5) - rscale * prop
	if robj < 0:
		return robj
	robj = math.pow(robj, 2.0/3.0)
	return robj

# Return the "time" coord in the Schwarzschild chart of a radially
# infalling massive object, currrently located at the radial
# coordinate `rcoord`.
def schw_tp_coord(rcoord):
	squ = math.sqrt(0.5*rcoord/mass)
	la = (squ - 1.0) / (squ + 1)
	if la < 0.0:
		la = -la
	lg = 2.0*mass* math.log(la)

	time = 2.0 * (rcoord + 6.0*mass) * squ / 3.0 + lg
	return time

# Return the "time" coord in the Schwarzschild chart of a radially
# upward-moving null geodesic, currrently located at the radial
# coordinate `rcoord`.
def schw_tnull_coord(rcoord):
	met = rcoord / (2 * mass) -1
	if met < 0:
		return rcoord - 2 * mass * math.log(-met)
	else:
		return rcoord + 2 * mass * math.log(met)

# Compute the difference in the time coordinate between the
# upward-moving null geodesic and the downward-moving astronaut.
# When this is zero, the two intersect.
# The goal is to solve for `rcoord` by using a root-finding algo.
# The second argument, rflash, is the location of the flashlight
# when it emits the signal.
def recv_null(rcoord, rflash):
	tastro = schw_tp_coord(rnaught) - schw_tp_coord(rcoord)
	tnull = schw_tnull_coord(rcoord) - schw_tnull_coord(rflash)
	diff = tnull - tastro
	return diff


tau = 0.0
tau = 38.0
taustep = 0.05
while True:

	# flashlight radial coordinate
	rflash = schw_rp_coord(rfnaught, tau)
	if rflash < 0.0:
		break

	# flashlight time coord
	tflash = schw_tp_coord(rfnaught) - schw_tp_coord(rflash)

	# astronaut radial coord, when astronaut receives the flash.
	rastro = scipy.optimize.brentq(recv_null, 0, 100, args=rflash)

	# astronaut time coord
	tastro = schw_tp_coord(rnaught) - schw_tp_coord(rastro)

	print( \
		"{:10.4f}".format(tau), \
		"{:10.4f}".format(rflash), \
		"{:10.4f}".format(tflash), \
		"{:10.4f}".format(rastro), \
		"{:10.4f}".format(tastro) \
		)
	tau += taustep
