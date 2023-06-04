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
# Note mass=1 corresponds to about 2e5 solar masses. 201900 to be exact.
mass = 1

# Initial location of lab frame
rnaught = 10

# Initial location of flashlight
deltar = -0.1

# Proper time step size
taustep = 0.5

tau = 0.0
#tau = 39.0     # when rnaught = 20
#taustep = 0.02

tau = 12.0      # when rnaught = 10
taustep = 0.01

# ------------------
# Solar mass, measured in seconds
mass = 5e-6
# Distance to Earth, measured in seconds
rnaught = 500

# Falling radially, from Earth to Sun, takes 27 days! (= 2.35e6 seconds)
tau = 0.0
tau = 2.35e6
taustep = 100
tau = 2356300.0
taustep = 0.1
tau = 2356315.0
taustep = 0.001
tau = 2356315.53
taustep = 1e-5

# Result: flashlight (located 0.1 in front) falls through EH at 2356315.5325
# and observer is still r=t=2.2 seconds away when this happens. (But this is
# a numerical artifact. See below.)

# If the infall is started at rnaught=1 and flashlight is 0.01 closer
# then flashlight falls in at 207.664145 while observer is still r=t=0.06 away.

# Curious: at this mass scale, there is both a separation in r and t
# from the EH, while, for the large masses, below, that r separation goes away.
# This hold true even if the flashlight is given a *huge* lead for the
# large masses.
#
# Answer: It is a numeric artifact. Both cross at the same r(=2m) but at
# different t coords. This is as it should be.

# ------------------
# Mass of Sgr A*, measured in seconds
mass = 20.5
# Distance to Earth, measured in seconds
rnaught = 500

tau = 0.0
taustep = 10
tau = 1110.0
taustep = 0.1
tau = 1135.0
taustep = 0.001

# Result: flashlight (located at 0.1 in front) falls through EH at 1136.3670
# Observer is r=0, t=0.19 seconds behind when this happens.
#
# But if flashlight is 1.0 in front, it falls in at 1133.2260
# Observer is r=0, t=1.9 seconds behind.
#
# But if flashlight is 10 in front, it falls in at 1101.97
# Observer is r=0, t=18.9 seconds behind.

# ------------------
# Mass of M87, measured in seconds
mass = 35000
rnaught = 100 * mass

tau = 0.0
taustep = 1000
tau = 16450000.0
taustep = 10
tau = 16452400.0
taustep = 1

# ------------------

mass = 1
rnaught = 50
deltar = -25

tau = 0.0
taustep = 1
tau = 57
taustep = 0.01
tau = 57.59
taustep = 1e-4
tau = 57.5922
taustep = 1e-6
tau = 57.592231
taustep = 1e-8
tau = 57.59223176
taustep = 1e-10
tau = 57.59223176
taustep = 1e-12

# ------------------
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
print("# astronaut r minus flashlight r")
print("# astronaut t minus flashlight t")
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
def recv_null(rcoord, rflash, tflash):
	tastro = schw_tp_coord(rnaught) - schw_tp_coord(rcoord)
	tnull = tflash + schw_tnull_coord(rcoord) - schw_tnull_coord(rflash)
	diff = tnull - tastro
	return diff

cross = False
while True:

	# flashlight radial coordinate
	rflash = schw_rp_coord(rfnaught, tau)
	if rflash < 0.0:
		break

	# Print a break when event horizon is crossed
	if cross == False and rflash < 2 * mass:
		cross = True
		print("# ----------------------------------------------------")
		break

	# flashlight time coord
	tflash = schw_tp_coord(rfnaught) - schw_tp_coord(rflash)

	# astronaut radial coord, when astronaut receives the flash.
	rmin = rflash
	rmax = rflash - 10*deltar
	rastro = scipy.optimize.brentq(recv_null, rmin, rmax, args=(rflash, tflash))

	# astronaut time coord
	tastro = schw_tp_coord(rnaught) - schw_tp_coord(rastro)

	print( \
		"{:10.12f}".format(tau), \
		"{:10.6f}".format(rflash), \
		"{:10.6f}".format(tflash), \
		"{:10.6f}".format(rastro), \
		"{:10.6f}".format(tastro), \
		"{:10.6f}".format(rastro-rflash), \
		"{:10.6f}".format(tastro-tflash) \
		)
	tau += taustep
