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
		fprintf(fpwrite, "Cuerpo: %d -> ", pos[i].id);
		fprintf(fpwrite, "\t%*f\t%*f\t%*f\t%*f\t%*f\t%*f\n", 10, pos[i].x, 10, pos[i].y, 10, vel[i].x, 10, vel[i].y, 10, acc[i].x, 10, acc[i].y);
	}

	fprintf(fpwrite, "\n");

}

void imprimirTerminal(int m){

	int i;
	for(i = 0; i < datos.n; i++){
		printf("Cuerpo: %d -> ", vel[i].id);
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

	int arrayNum[4] = {1, 2, 3, 4};
	int *array;
	array = arrayNum;

	if(rank == 0){
	 	
		
	// 	// (array[0])=0;
	// 	// (array[1])=1;
	// 	// (array[2])=2;
	// 	// (array[3])=3;
	}else{
	 	
	}
	array = malloc(sizeof(int)*npr*2/npr);

	MPI_Scatter(arrayNum, 4/npr, MPI_INT, array, 4/npr, MPI_INT, 0, MPI_COMM_WORLD);

	for(int i=0; i<2; i++) {
		printf("Rank %d: a[%d]=%d\n", rank, i, *(array + i));
	}

	MPI_Finalize();


	return 0;

}
