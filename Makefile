


FC=gfortran

all: 0llin

0llin: src/main.f90
	$(FC) src/main.f90 -o 0llin

