#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <openmpi/mpi.h>
#include "timer.h"

#define FREAD "datos.dat"
#define FWRITE "out.txt"
#define G 1

FILE *fpread, *fpwrite;

struct Datos {
	int n;
	int tp;
	int k;
	double delta;
	double u;
};

struct Masas {
	int id;
	double m;
};

struct Coord {
	int id;
	double x;
	double y;
};

struct Datos datos;
struct Masas *masas;
struct Coord *pos; // Posiciones
struct Coord *vel; // Velocidades
struct Coord *acc; // Aceleraciones

// MPI Types
MPI_Datatype MPI_Datos;
MPI_Datatype MPI_Masas;
MPI_Datatype MPI_Coord;

// Variable cuerpos totales
int cuerpos_totales;

void imprimirFichero(){

	int i;
	for(i = 0; i < datos.n; i++){
		fprintf(fpwrite, "Cuerpo: %d -> ", i);
		fprintf(fpwrite, "\t%*f\t%*f\t%*f\t%*f\t%*f\t%*f\n", 10, pos[i].x, 10, pos[i].y, 10, vel[i].x, 10, vel[i].y, 10, acc[i].x, 10, acc[i].y);
	}

	fprintf(fpwrite, "\n");

}

void imprimirTerminal(int m){

	int i;
	for(i = 0; i < datos.n; i++){
		printf("Cuerpo: %d -> ", i);
		if(m) printf("Masa: %.2f\n", masas[i].m);
		printf("\t%*f\t%*f\t%*f\t%*f\t%*f\t%*f\n", 10, pos[i].x, 10, pos[i].y, 10, vel[i].x, 10, vel[i].y, 10, acc[i].x, 10, acc[i].y);
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
		fscanf(fpread, "%lf, %lf, %lf, %lf, %lf", &(masas[i]).m, &(pos[i]).x, &(pos[i]).y, &(vel[i]).x, &(vel[i]).y);
		masas[i].id = i; pos[i].id = i; vel[i].id = i;
		acc[i].x = 0; acc[i].y = 0;
	}
	// Cuerpos vacíos
	for(i = datos.n; i < cuerpos_totales; i++){
		masas[i].id = -1; masas[i].m = 0;
		pos[i].id = -1; pos[i].x = 0; pos[i].y = 0;
		vel[i].id = -1; vel[i].x = 0; vel[i].y = 0;
		acc[i].id = -1; acc[i].x = 0; acc[i].y = 0;
	}

	fclose(fpread);

}

void calcularAceleracion(){/*

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

*/}

void mpi_datatype_datos(){

	int blcklen[2] = {3, 2};
	MPI_Aint displ[2] = {offsetof(struct Datos, n), offsetof(struct Datos, delta)};
	MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
	MPI_Type_create_struct(2, blcklen, displ, types, &MPI_Datos);
	MPI_Type_commit(&MPI_Datos);

}

void mpi_datatype_masas(){

	int blcklen[2] = {1, 1};
	MPI_Aint displ[2] = {offsetof(struct Masas, id), offsetof(struct Masas, m)};
	MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
	MPI_Type_create_struct(2, blcklen, displ, types, &MPI_Masas);
	MPI_Type_commit(&MPI_Masas);

}

void mpi_datatype_coord(){

	int blcklen[2] = {1, 2};
	MPI_Aint displ[2] = {offsetof(struct Coord, id), offsetof(struct Coord, x)};
	MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
	MPI_Type_create_struct(2, blcklen, displ, types, &MPI_Coord);
	MPI_Type_commit(&MPI_Coord);

}

int main(int argc, char *argv[]) {

	/* - - - - - Inicialización - - - - - */

	int rank, npr;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &npr);

	// Creación datatypes MPI

	mpi_datatype_datos();
	mpi_datatype_masas();
	mpi_datatype_coord();

	// Rank 0 obtiene datos de entrada

	if(rank == 0){

			int opcion;

			// PREGUNTA: Lectura de datos.

			printf("¿Cómo debe obtener el programa los datos de entrada?\n (1) A través del fichero \"datos.dat\".\n (2) Por teclado.\n");
			scanf("%d", &opcion);

			while(opcion != 1 && opcion != 2){
				printf("ERROR. Opción no disponible.\n Vuelva a introducir el dato.\n (1) A través del fichero \"datos.dat\".\n (2) Por teclado.\n");
				scanf("%d", &opcion);
			}

			if(opcion == 1) leerFichero();
			else leerTeclado();

	}

	// Envío/Recepción de la estructura con los datos

	MPI_Bcast(&datos, 1, MPI_Datos, 0, MPI_COMM_WORLD);

	// Se calculan los cuerpos totales (contando vacíos)

	int ncu = datos.n / npr;
	if((datos.n % npr) > 0) ncu++;
	cuerpos_totales = ncu * npr;

	// Se reserva memoria para guardar los datos de los cuerpos

	int aux = (rank == 0)?cuerpos_totales:ncu;
	masas = malloc(sizeof(struct Masas) * aux);
	pos = malloc(sizeof(struct Coord) * cuerpos_totales);
	vel = malloc(sizeof(struct Coord) * aux);
	acc = malloc(sizeof(struct Coord) * aux);

	// Rank 0 obtiene y muestra/escribe los datos de los cuerpos

	if(rank == 0){

			int opcion;

			// PREGUNTA: Mostrar salida.

			printf("\n¿Cómo desea visualizar los datos resultantes del programa?\n (1) A través de la terminal.\n (2) En un fichero de texto.\n");
			scanf("%d", &opcion);

			while(opcion != 1 && opcion != 2){
				printf("ERROR. Opción no disponible.\n Vuelva a introducir el dato.\n (1) A través de la terminal.\n (2) En un fichero de texto.\n");
				scanf("%d", &opcion);
			}

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

	}

	// Envío/Recepción de las posiciones de los cuerpos

	MPI_Scatter(vel, ncu, MPI_Coord, vel, ncu, MPI_Coord, 0, MPI_COMM_WORLD);

	// if(rank==0) for(int i = 0; i < ncu; i++) printf("%d: %f - %f\n", rank, vel[i].x, vel[i].y); // DEBUG
	// if(rank==1) for(int i = 0; i < ncu; i++) printf("%d: %f - %f\n", rank, vel[i].x, vel[i].y); // DEBUG

	// Envío/Recepción de las posiciones de los cuerpos

	MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, pos, cuerpos_totales, MPI_Coord, MPI_COMM_WORLD);

	// if(rank == 0) for(int i = 0; i < cuerpos_totales; i++) printf("%d: %f - %f\n", rank, pos[i].x, pos[i].y); // DEBUG
	// if(rank == 1) for(int i = 0; i < cuerpos_totales; i++) printf("%d: %f - %f\n", rank, pos[i].x, pos[i].y); // DEBUG

	MPI_Finalize();

	free(masas);
	free(pos);
	free(vel);
	free(acc);

	return 0;

}
