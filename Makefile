all: compile-par run-par

compile-sec:
	gcc NCuerposSecuencial.c -o NCuerposSecuencial -lm -Wall

run-sec:
	./NCuerposSecuencial

compile-par:
	mpicc NCuerposParalelo.c -o NCuerposParalelo -lm -Wall

compile-par-nosal:
	mpicc NCuerposParalelo.c -o NCuerposParalelo -lm -Wall -D NO_SAL

run-par:
	mpirun -np 2 NCuerposParalelo
