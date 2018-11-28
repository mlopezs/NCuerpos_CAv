#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "timer.h"

#define FICHERO "datos.dat"

struct DatosCuerpo {
	double masa;
	double posicionX;
	double posicionY;
	double velocidadX;
	double velocidadY;
	double aceleracionX;
	double aceleracionY;
};

FILE *fichero;
int n, tp, k;
double delta, u;
struct DatosCuerpo *cuerpos;

void leerFichero()
{
	
	/*---------Lectura del archivo--------*/

		if((fichero = fopen(FICHERO,"r")) == NULL) {
			printf("Error al abrir el archivo");
			exit(EXIT_FAILURE);
		}

		fscanf(fichero, "%d, %lf, %d, %lf, %d", &n, &delta, &tp, &u, &k);

		printf("N: %d\n", n);
		printf("Delta: %.2f\n", delta);
		printf("Tp: %d\n", tp);
		printf("U: %.2f\n", u);
		printf("K: %d\n", k);
		
		fclose(fichero);
}

void leerEntradaTeclado()
{
	printf("Introduzca el número de cuerpos (n) para el programa\n");
	scanf("%d", &n);
	printf("Introduzca el incremento del tiempo en cada paso (delta)\n");
	scanf("%lff", &delta);
	printf("Introduzca el número total de pasos (tp)\n");
	scanf("%d", &tp);
	printf("Introduzca la distancia umbral (u)\n");
	scanf("%lff", &u);
	printf("Introduzca el k deseado\n");
	scanf("%d", &k);

	printf("Los datos introducidos son: %d, %.2f, %d, %.2f, %d\n", n, delta, tp, u, k);
}

void leerDatosCuerpo()
{
	int i;
	if((fichero = fopen(FICHERO,"r")) == NULL) {
		printf("Error al abrir el archivo");
		exit(EXIT_FAILURE);
	}

	fscanf(fichero, "%*[^\n]\n"); /*--- Para situarnos en la segunda línea del archivo ---*/
	for (i = 0; i < n; i++)
	{
		fscanf(fichero, "%lf, %lf, %lf, %lf, %lf", &(cuerpos[i]).masa, &(cuerpos[i]).posicionX, &(cuerpos[i]).posicionY,
			&(cuerpos[i]).aceleracionX, &(cuerpos[i]).aceleracionY);

		printf("Los datos leídos son: Masa: %.1f " "PosicionX: %.3f " "PosicionY: %.3f " "AceleracionX: %.3f " "AceleracionY: %.3f\n",
			cuerpos[i].masa, cuerpos[i].posicionX, cuerpos[i].posicionY, cuerpos[i].aceleracionX, cuerpos[i].aceleracionY);
	}
	fclose(fichero);
}


int main(int argc, char const *argv[])
{
	int opcionElegida;
	double inicio, fin;
	
	printf("¿Cómo desea leer los datos para el programa?\n -1. A través de fichero.\n -2. Por entrada de teclado.\n");
	scanf("%d", &opcionElegida);

	if (opcionElegida == 1)
	{
		leerFichero();	
	} else {
		leerEntradaTeclado();
	}

	/*---------- Reserva de espacio para los n cuerpos ---------*/
	cuerpos = malloc(sizeof(struct DatosCuerpo)*n);
	
	/*---------- Lectura de masa, posiciones y velocidades iniciales ------*/
	leerDatosCuerpo();


	return 0;
}
