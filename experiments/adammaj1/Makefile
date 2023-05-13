
all: tav elec

elec: electric.c
	cc -O2 $< -o elec -lm

tav: tavis.cpp
	cc -O2 $< -o tav -lm -lstdc++
