#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "timer.h"

#define FICHERO "datos.dat"

#define G 1

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

		printf("Los datos leídos del fichero son: %d, %.2f, %d, %.2f, %d\n", n, delta, tp, u, k);

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
	char buffer[1024];
	if((fichero = fopen(FICHERO,"r")) == NULL) {
		printf("Error al abrir el archivo");
		exit(EXIT_FAILURE);
	}

	fgets(buffer, 1024, fichero); /*--- Para situarnos en la segunda línea del archivo ---*/
	for (i = 0; i < n; i++)
	{
		fscanf(fichero, "%lf, %lf, %lf, %lf, %lf", &(cuerpos[i]).masa, &(cuerpos[i]).posicionX, &(cuerpos[i]).posicionY,
			&(cuerpos[i]).velocidadX, &(cuerpos[i]).velocidadY);

		printf("Los datos leídos son: Masa: %.1f " "PosicionX: %.3f " "PosicionY: %.3f " "VelocidadX: %.3f " "VelocidadY: %.3f\n",
			cuerpos[i].masa, cuerpos[i].posicionX, cuerpos[i].posicionY, cuerpos[i].velocidadX, cuerpos[i].velocidadY);
	}
	fclose(fichero);
}

void calcularAceleracion(){

	int i;
	for(i = 0; i < n; i++){
		cuerpos[i].aceleracionX = 0.0;
		cuerpos[i].aceleracionY = 0.0;
	}

	int q, p;
	double distX, distY, distMod;
	double dist3;
	double accX, accY;

	for(q = 0; q < n; q++){
		for(p = q+1; p < n; p++){

			distX = cuerpos[q].posicionX - cuerpos[p].posicionX;
			distY = cuerpos[q].posicionY - cuerpos[p].posicionY;
			distMod = sqrt(pow(distX,2) + pow(distY,2));

			if(distMod < u){
				printf("Distancia menor que el umbral.\n");
			}else{

				dist3 = pow(distMod, 3);

				distX = cuerpos[p].posicionX - cuerpos[q].posicionX;
				distY = cuerpos[p].posicionY - cuerpos[q].posicionY;

				accX = (G * distX) / dist3;
				accY = (G * distY) / dist3;

				cuerpos[q].aceleracionX += accX * cuerpos[p].masa;
				cuerpos[q].aceleracionY += accY * cuerpos[p].masa;

				cuerpos[p].aceleracionX += accX * -cuerpos[q].masa;
				cuerpos[p].aceleracionY += accY * -cuerpos[q].masa;

			}
		}
	}

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

	/*---------- Algoritmo ----------*/

	int i, paso, q;

	GET_TIME(inicio);

	calcularAceleracion();

	printf("Aceleraciones iniciales:\n");
	for(i = 0; i < n; i++){
		printf("AceleracionX: %.3f, AceleracionY: %.3f\n", cuerpos[i].aceleracionX, cuerpos[i].aceleracionY);
	}

	int t = 0;

	for(paso = 1; paso <= tp; paso++){
		for(q = 0; q < n; q++){
			cuerpos[q].posicionX += cuerpos[q].velocidadX * delta;
			cuerpos[q].posicionY += cuerpos[q].velocidadY * delta;
			cuerpos[q].velocidadX += cuerpos[q].aceleracionX * delta;
			cuerpos[q].velocidadY += cuerpos[q].aceleracionY * delta;

			printf("PosicionX: %.3f, PosicionY: %.3f, VelocidadX: %.3f, VelocidadY: %.3f\n", cuerpos[q].posicionX, cuerpos[q].posicionY, cuerpos[q].velocidadX, cuerpos[q].velocidadY);
		}
		for(q = 0; q < n; q++){
			calcularAceleracion();
			printf("Aceleraciones:\n");
			printf("AceleracionX: %.3f, AceleracionY: %.3f\n", cuerpos[q].aceleracionX, cuerpos[q].aceleracionY);
		}
		t += delta;
		printf("%.2d\n", t);
	}

	GET_TIME(fin);
	printf("\nEjecución en %.2f segundos.\n", fin - inicio);

	free(cuerpos);

	return 0;
}
