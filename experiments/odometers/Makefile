
EXES= bernoulli-odo ornstein-odo

all: ${EXES}

bernoulli-odo: bernoulli-odo.o
ornstein-odo: ornstein-odo.o

${EXES}:
	cc -o $@ $< -lm

.c.o:
	cc -c -O2 $<
