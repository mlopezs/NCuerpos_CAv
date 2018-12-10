all: compile run

compile:
	gcc NCuerposSecuencial.c -o NCuerposSecuencial -lm -Wall

run:
	./NCuerposSecuencial 
