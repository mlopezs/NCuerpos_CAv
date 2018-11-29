#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "timer.h"

#define FICHERO "datos.dat"
#define G 1

struct DatosCuerpo {
	int id;
	double masa;
	double posicionX;
	double posicionY;
	double velocidadX;
	double velocidadY;
	double aceleracionX;
	double aceleracionY;
};

FILE *fichero, *fp;
int n, tp, k;
double delta, u;
struct DatosCuerpo *cuerpos;

void imprimirCuerpos(int m, FILE *fp){
	int i;
	
	for(i = 0; i < n; i++){
		fprintf(fp, "Cuerpo: %d -> ", cuerpos[i].id);
		if(m) printf("Masa: %.2f\n", cuerpos[i].masa);
		/*fprintf(fp, "Pos: (%.10f,%10f), Vel: (%.10f,%.10f), Ace: (%.10f,%.10f)\n", cuerpos[i].posicionX, 
			cuerpos[i].posicionY, cuerpos[i].velocidadX, cuerpos[i].velocidadY, cuerpos[i].aceleracionX, cuerpos[i].aceleracionY);*/

		fprintf(fp, "  %f     %f     %f      %f      %f        %f\n", cuerpos[i].posicionX, cuerpos[i].posicionY, cuerpos[i].velocidadX, 
				cuerpos[i].velocidadY, cuerpos[i].aceleracionX, cuerpos[i].aceleracionY);
	}
	fprintf(fp, "\n");
}

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
		cuerpos[i].id = i;

		printf("Los datos leídos son: Masa: %.1f " "PosicionX: %.3f " "PosicionY: %.3f " "VelocidadX: %.3f " "VelocidadY: %.3f\n",
			cuerpos[i].masa, cuerpos[i].posicionX, cuerpos[i].posicionY, cuerpos[i].velocidadX, cuerpos[i].velocidadY);
	}

	imprimirCuerpos(1, fp);

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
				/*printf("Distancia menor que el umbral.\n");*/
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


	fp = fopen( "out_file.txt", "w" ); /*--- Open file for writing */
	fprintf(fp, "Por cada instante de tiempo y cada cuerpo aparecen: \n");
	fprintf(fp, "\t      Posicion(x), Posicion(y), Velocidad(x), Velocidad(y), Aceleracion(x), Aceleracion(y)\n \n");
	fprintf(fp, "0.00\n");

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

	int paso, q;

	GET_TIME(inicio);

	double t = 0.0;
	
	calcularAceleracion();

	for(paso = 1; paso <= tp; paso++){
		for(q = 0; q < n; q++){
			cuerpos[q].posicionX += cuerpos[q].velocidadX * delta;
			cuerpos[q].posicionY += cuerpos[q].velocidadY * delta;
			cuerpos[q].velocidadX += cuerpos[q].aceleracionX * delta;
			cuerpos[q].velocidadY += cuerpos[q].aceleracionY * delta;
		}
		for(q = 0; q < n; q++){
			calcularAceleracion();
		}

		t += delta;

		if(paso % k == 0){
			fprintf(fp, "%.2f\n", t);
			imprimirCuerpos(0, fp);
			printf("\n");
		}


	}

	GET_TIME(fin);
	
	printf("\nEjecución en %f segundos.\n", (fin - inicio));

	free(cuerpos);

	return 0;
}
