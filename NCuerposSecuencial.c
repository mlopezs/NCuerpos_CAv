#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "timer.h"

#define FREAD "datos.dat"
#define FWRITE "out.txt"
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

FILE *fpread, *fpwrite;
int n, tp, k;
double delta, u;
struct DatosCuerpo *cuerpos;

void imprimirFichero(){

	int i;
	for(i = 0; i < n; i++){
		fprintf(fpwrite, "Cuerpo: %d -> ", cuerpos[i].id);
		fprintf(fpwrite, "\t%*f\t%*f\t%*f\t%*f\t%*f\t%*f\n", 10, cuerpos[i].posicionX, 10, cuerpos[i].posicionY, 10, cuerpos[i].velocidadX,
				10, cuerpos[i].velocidadY, 10, cuerpos[i].aceleracionX, 10, cuerpos[i].aceleracionY);
	}

	fprintf(fpwrite, "\n");

}

void imprimirTerminal(int m){
	int i;
	for(i = 0; i < n; i++){
		printf("Cuerpo: %d -> ", cuerpos[i].id);
		if(m) printf("Masa: %.2f\n", cuerpos[i].masa);
		printf("\t%*f\t%*f\t%*f\t%*f\t%*f\t%*f\n", 10, cuerpos[i].posicionX, 10, cuerpos[i].posicionY, 10, cuerpos[i].velocidadX, 10, cuerpos[i].velocidadY, 10, cuerpos[i].aceleracionX, 10, cuerpos[i].aceleracionY);
	}
}

void leerFichero() {

	if((fpread = fopen(FREAD,"r")) == NULL) {
		printf("Error al abrir el archivo");
		exit(EXIT_FAILURE);
	}

	fscanf(fpread, "%d, %lf, %d, %lf, %d", &n, &delta, &tp, &u, &k);

	printf("Los datos leídos del fichero son: n=%d, delta=%.2f, tp=%d, u=%.2f, k=%d\n", n, delta, tp, u, k);

	fclose(fpread);
}

void leerTeclado() {

	printf("Introduzca el número de cuerpos (n) para el programa: \n");
	scanf("%d", &n);
	printf("Introduzca el incremento del tiempo (delta) en cada paso: \n");
	scanf("%lff", &delta);
	printf("Introduzca el número total de pasos (tp): \n");
	scanf("%d", &tp);
	printf("Introduzca la distancia umbral (u): \n");
	scanf("%lff", &u);
	printf("Introduzca el k deseado: \n");
	scanf("%d", &k);

}

void leerDatosCuerpo() {

	if((fpread = fopen(FREAD,"r")) == NULL) {
		printf("ERROR. La apertura del archivo falló.");
		exit(EXIT_FAILURE);
	}

	// Nos colocamos en la segunda línea
	char buffer[1024];
	fgets(buffer, 1024, fpread);

	int i;
	for (i = 0; i < n; i++) {
		fscanf(fpread, "%lf, %lf, %lf, %lf, %lf", &(cuerpos[i]).masa, &(cuerpos[i]).posicionX, &(cuerpos[i]).posicionY,
		&(cuerpos[i]).velocidadX, &(cuerpos[i]).velocidadY);
		cuerpos[i].id = i;
	}

	fclose(fpread);

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

			if(distMod >= u){ // Control umbral

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

int main(int argc, char const *argv[]) {

	int opcion;

	/* - - - - - PREGUNTA: Lectura de datos. - - - - - */

	printf("¿Cómo debe obtener el programa los datos de entrada?\n (1) A través del fichero \"datos.dat\".\n (2) Por teclado.\n");
	scanf("%d", &opcion);

	while(opcion != 1 && opcion != 2){
		printf("ERROR. Opción no disponible.\n Vuelva a introducir el dato.\n (1) A través del fichero \"datos.dat\".\n (2) Por teclado.\n");
		scanf("%d", &opcion);
	}

	if(opcion == 1) leerFichero();
	else leerTeclado();

	/* - - - - - PREGUNTA: Mostrar salida. - - - - - */

	printf("\n¿Cómo desea visualizar los datos resultantes del programa?\n (1) A través de la terminal.\n (2) En un fichero de texto.\n");
	scanf("%d", &opcion);

	while(opcion != 1 && opcion != 2){
		printf("ERROR. Opción no disponible.\n Vuelva a introducir el dato.\n (1) A través de la terminal.\n (2) En un fichero de texto.\n");
		scanf("%d", &opcion);
	}

	/* - - - - - Comienzo del programa - - - - - */

	if(opcion == 2){
		fpwrite = fopen( FWRITE, "w" );
		fprintf(fpwrite, "Por cada instante de tiempo y cada cuerpo aparecen: \n");
		fprintf(fpwrite, "            \t%*s \t%*s \t%*s \t%*s \t%*s \t%*s\n", 10, "Posicion(x)", 10, "Posicion(y)", 10, "Velocidad(x)", 10, "Velocidad(y)", 10, "Aceleracion(x)", 10, "Aceleracion(y)");
	}

	// Reserva de memoria para 'n' cuerpos
	cuerpos = malloc( sizeof(struct DatosCuerpo) * n);

	// Lectura de datos de cada cuerpo (Siempre por fichero)
	leerDatosCuerpo();

	// Muestra de los datos
	printf("\nDescripción de los cuerpos:\n");
	imprimirTerminal(1);

	/* - - - - - Comienzo del algoritmo - - - - - */

	int paso, q;
	double inicio, fin;
	int flag = 1;

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
			if(opcion == 1){
				if(flag){
					 printf("\n            \t%*s\t%*s\t%*s\t%*s\t%*s\t%*s\n", 10, "Posicion(x)", 10, "Posicion(y)", 10, "Velocidad(x)", 10, "Velocidad(y)", 10, "Aceleracion(x)", 10, "Aceleracion(y)");
					 flag = 0;
				 }
	 			printf("%.2f\n", t);
				imprimirTerminal(0);
			}else{
				fprintf(fpwrite, "%.2f\n", t);
				imprimirFichero();
			}
		}
	}

	GET_TIME(fin);

	printf("\nEjecución en %f segundos.\n", (fin - inicio));
	if(opcion == 2) fprintf(fpwrite, "Programa ejecutado en %f segundos.\n", (fin - inicio));

	free(cuerpos);

	return 0;
}
