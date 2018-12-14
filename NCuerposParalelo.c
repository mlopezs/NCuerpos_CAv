#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <openmpi/mpi.h>
#include "timer.h"

#define FREAD "datos.dat"
#define FWRITE "out.txt"
#define G 1

struct Datos {
	int n;
	int tp;
	int k;
	double delta;
	double u;
};

struct Cuerpo {
	int id;
	double masa;
	double posX;
	double posY;
	double velX;
	double velY;
	double accX;
	double accY;
};

struct Datos datos;
struct Cuerpo *cuerpos;

FILE *fpread, *fpwrite;

void imprimirFichero(){

	int i;
	for(i = 0; i < datos.n; i++){
		fprintf(fpwrite, "Cuerpo: %d -> ", cuerpos[i].id);
		fprintf(fpwrite, "\t%*f\t%*f\t%*f\t%*f\t%*f\t%*f\n", 10, cuerpos[i].posX, 10, cuerpos[i].posY, 10, cuerpos[i].velX,
				10, cuerpos[i].velY, 10, cuerpos[i].accX, 10, cuerpos[i].accY);
	}

	fprintf(fpwrite, "\n");

}

void imprimirTerminal(int m){

	int i;
	for(i = 0; i < datos.n; i++){
		printf("Cuerpo: %d -> ", cuerpos[i].id);
		if(m) printf("Masa: %.2f\n", cuerpos[i].masa);
		printf("\t%*f\t%*f\t%*f\t%*f\t%*f\t%*f\n", 10, cuerpos[i].posX, 10, cuerpos[i].posY, 10, cuerpos[i].velX, 10, cuerpos[i].velY, 10, cuerpos[i].accX, 10, cuerpos[i].accY);
	}

}

void leerFichero() {

	if((fpread = fopen(FREAD,"r")) == NULL) {
		printf("Error al abrir el archivo");
		exit(EXIT_FAILURE);
	}

	fscanf(fpread, "%d, %lf, %d, %lf, %d", &datos.n, &datos.delta, &datos.tp, &datos.u, &datos.k);

	printf("Los datos leídos del fichero son: n=%d, delta=%.2f, tp=%d, u=%.2f, k=%d\n", datos.n, datos.delta, datos.tp, datos.u, datos.k);

	fclose(fpread);

}

void leerTeclado() {

	printf("Introduzca el número de cuerpos (n) para el programa: \n");
	scanf("%d", &datos.n);
	printf("Introduzca el incremento del tiempo (delta) en cada paso: \n");
	scanf("%lff", &datos.delta);
	printf("Introduzca el número total de pasos (tp): \n");
	scanf("%d", &datos.tp);
	printf("Introduzca la distancia umbral (u): \n");
	scanf("%lff", &datos.u);
	printf("Introduzca el k deseado: \n");
	scanf("%d", &datos.k);

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
	for (i = 0; i < datos.n; i++) {
		fscanf(fpread, "%lf, %lf, %lf, %lf, %lf", &(cuerpos[i]).masa, &(cuerpos[i]).posX, &(cuerpos[i]).posY, &(cuerpos[i]).velX, &(cuerpos[i]).velY);
		cuerpos[i].id = i;
	}

	fclose(fpread);

}

void calcularAceleracion(){

	int i;
	for(i = 0; i < datos.n; i++){
		cuerpos[i].accX = 0.0;
		cuerpos[i].accY = 0.0;
	}

	int q, p;
	double distX, distY, distMod;
	double dist3;
	double accX, accY;

	for(q = 0; q < datos.n; q++){
		for(p = q+1; p < datos.n; p++){

			distX = cuerpos[q].posX - cuerpos[p].posX;
			distY = cuerpos[q].posY - cuerpos[p].posY;
			distMod = sqrt(pow(distX,2) + pow(distY,2));

			if(distMod >= datos.u){ // Control umbral

				dist3 = pow(distMod, 3);

				distX = cuerpos[p].posX - cuerpos[q].posX;
				distY = cuerpos[p].posY - cuerpos[q].posY;

				accX = (G * distX) / dist3;
				accY = (G * distY) / dist3;

				cuerpos[q].accX += accX * cuerpos[p].masa;
				cuerpos[q].accY += accY * cuerpos[p].masa;

				cuerpos[p].accX += accX * -cuerpos[q].masa;
				cuerpos[p].accY += accY * -cuerpos[q].masa;

			}
		}
	}

}

int main(int argc, char *argv[]) {

	/* - - - - - Inicialización MPI - - - - - */

	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// Creamos un MPI datatype para los datos
	MPI_Datatype MPI_Datos;
	int blcklen[2] = {3, 2};
	MPI_Aint displ[2] = {offsetof(struct Datos, n), offsetof(struct Datos, delta)};
	MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
	MPI_Type_create_struct(2, blcklen, displ, types, &MPI_Datos);
	MPI_Type_commit(&MPI_Datos);

	if(rank == 0){ // Proceso con rango 0

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

	}

	// Reserva de memoria para 'n' cuerpos
	cuerpos = malloc( sizeof(struct Cuerpo) * datos.n);

	if(rank == 0){

			int opcion;

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

			// Lectura de datos de cada cuerpo (Siempre por fichero)
			leerDatosCuerpo();

			// Muestra de los datos
			printf("\nDescripción de los cuerpos:\n");
			imprimirTerminal(1);


			// Enviar datos a los demas procesos
			MPI_Bcast(&datos, 1, MPI_Datos, 0, MPI_COMM_WORLD);

			// Enviar los datos de los cuerpos a los demas procesos


			/* - - - - - Comienzo del algoritmo - - - - - */

			int paso, q;
			double inicio, fin;
			int flag = 1;

			GET_TIME(inicio);

			double t = 0.0;

			calcularAceleracion();

			for(paso = 1; paso <= datos.tp; paso++){

				for(q = 0; q < datos.n; q++){
					cuerpos[q].posX += cuerpos[q].velX * datos.delta;
					cuerpos[q].posY += cuerpos[q].velY * datos.delta;
					cuerpos[q].velX += cuerpos[q].accX * datos.delta;
					cuerpos[q].velY += cuerpos[q].accY * datos.delta;
				}

				for(q = 0; q < datos.n; q++){
					calcularAceleracion();
				}

				t += datos.delta;

				if(paso % datos.k == 0){
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

			// Fin main

	} else { // Procesos con rango > 0

		// Reciben los datos del proceso 0
		MPI_Bcast(&datos, 1, MPI_Datos, 0, MPI_COMM_WORLD);

	}

	// Creamos un MPI datatype para los datos de cuerpos
	MPI_Datatype MPI_Cuerpos;
	int blcklen2[2] = {1, 7};
	MPI_Aint displ2[2] = {offsetof(struct Cuerpo, id), offsetof(struct Cuerpo, posX)};
	MPI_Datatype types2[2] = {MPI_INT, MPI_DOUBLE};
	MPI_Type_create_struct(datos.n, blcklen2, displ2, types2, &MPI_Cuerpos);
	MPI_Type_commit(&MPI_Cuerpos);

	if(rank == 0){

		MPI_Bcast(&cuerpos, 1, MPI_Cuerpos, 0, MPI_COMM_WORLD);

	} else {

		// Recibe los datos de los cuerpos del proceso 0
		MPI_Bcast(&cuerpos, 1, MPI_Cuerpos, 0, MPI_COMM_WORLD);

		/*




		Comprobar si se envia y recibe bien el cuerpos







		*/

	}

	MPI_Finalize();

	return 0;
}
