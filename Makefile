all: compile-par-rapido run-par-rapido

compile-sec:
	mpicc NCuerposSecuencial.c -o NCuerposSecuencial -lm -Wall

run-sec:
	./NCuerposSecuencial

compile-par:
	mpicc NCuerposParalelo.c -o NCuerposParalelo -lm -Wall

compile-par-nosal:
	mpicc NCuerposParalelo.c -o NCuerposParalelo -lm -Wall -D NO_SAL

run-par:
	mpirun -np 2 NCuerposParalelo

compile-par-rapido:
	mpicc NCuerposParalelo_AlgoritmoRapido.c -o NCuerposParalelo_AlgoritmoRapido -lm -Wall

compile-par-rapido-nosal:
	mpicc NCuerposParalelo_AlgoritmoRapido.c -o NCuerposParalelo_AlgoritmoRapido -lm -Wall -D NO_SAL

run-par-rapido:
	mpirun -np 2 NCuerposParalelo_AlgoritmoRapido
