all: compile-par run-par

compile-sec:
	gcc NCuerposSecuencial.c -o NCuerposSecuencial -lm -Wall

run-sec:
	./NCuerposSecuencial

compile-par:
	mpicc NCuerposParalelo.c -o NCuerposParalelo -lm -Wall

run-par:
	mpirun -np 2 NCuerposParalelo
