#! /usr/bin/env python3
#
# Solve the equations for a radially infalling flashlight,
# followed by a radially infalling astronaut. Record what
# the coordinates of each. All values are relative to the
# proper time of the flashlight, when the flashlight flashed.
#

mass = 1
rnaught = 20
deltar = -0.1

taustep = 0.1

print("#\n# Radial infall simulation")
print("# Mass=", mass)
print("# r_0=", rnaught)
print("# delta r_0=", deltar)

tau = 0.0
while tau < 10.0:
	print("duuude", "{:10.4f}".format(tau))
	tau += taustep
