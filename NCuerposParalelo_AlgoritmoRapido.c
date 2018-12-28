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
struct Coord *pos; // Posiciones de MIS cuerpos
struct Coord *vel; // Velocidades
struct Coord *acc; // Aceleraciones

struct Coord *p_local;
struct Coord *a_local;
struct Coord *p_anillo;
struct Coord *a_anillo;

// MPI Types
MPI_Datatype MPI_Datos;
MPI_Datatype MPI_Masas;
MPI_Datatype MPI_Coord;

MPI_Datatype MPI_CNCR;

int rank, npr; // Id nodo y numero de procesadores

int cuerpos_totales; // Cuerpos totales (+vacíos)
int ncu; // Cuerpos por procesador

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

void calcularAceleracion(int fase){

	double dx, dy, dm, d3, ax, ay;

	int yo = (rank + fase) % npr;
	int rk = rank;

	for(int c = 0; c < ncu; c++){
		for(int d = 0; d < ncu; d++){

			if(yo > rk){

				dx = p_local[c].x - p_local[d].x;
				dy = p_local[c].y - p_local[d].y;
				dm = sqrt(pow(dx, 2) + pow(dy, 2));

				if(dm >= datos.u){

					d3 = pow(dm, 3);

					ax = (G * dx) / d3;
					ay = (G * dy) / d3;

					a_local[d].x += ax * masas[p_local[c].id].m;
					a_local[d].y += ay * masas[p_local[c].id].m;

					a_anillo[c].x += ax * -masas[p_local[d].id].m;
					a_anillo[c].y += ay * -masas[p_local[d].id].m;

				}
			}
			rk += npr;
		}
		yo += npr;
	}

}

void prepararAceleracion(){

	int src = (rank + 1) % npr;
	int dst = (rank - 1 + npr) % npr;

	// p_anillo <- p_local

	for(int i = 0; i < ncu; i++){
		p_anillo[i].id = p_local[i].id;
		p_anillo[i].x = p_local[i].x;
		p_anillo[i].y = p_local[i].y;
	}

	// a_local <- a_anillo <- 0.0

	for(int i = 0; i < ncu; i++){
		a_anillo[i].id = p_local[i].id;
		a_anillo[i].x = 0.0;
		a_anillo[i].y = 0.0;
		a_local[i].id = p_local[i].id;
		a_local[i].x = 0.0;
		a_local[i].y = 0.0;
	}

	// for(int i = 0; i < ncu; i++){
	// 	a_anillo[i].id = a_local[i].id = p_local[i].id;
	// 	a_local[i].x = a_anillo[i].x = 0.0;
	// 	a_local[i].y = a_anillo[i].y = 0.0;
	// }

	// Calculamos aceleraciones

	calcularAceleracion(0);

	for(int fase = 1; fase < npr; fase++){
		MPI_Sendrecv_replace(p_anillo, ncu, MPI_Coord, dst, 1, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		MPI_Sendrecv_replace(a_anillo, ncu, MPI_Coord, dst, 1, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		calcularAceleracion(fase);

	}

	MPI_Sendrecv_replace(a_anillo, ncu, MPI_Coord, dst, 1, src, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	for(int i = 0; i < ncu; i++){
		a_local[i].x += a_anillo[i].x;
		a_local[i].y += a_anillo[i].y;
	}

}



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

void mpi_datatype_cncr(){
	MPI_Datatype MPI_CNC;
	MPI_Type_vector(ncu, 1, npr, MPI_Coord, &MPI_CNC);
	MPI_Type_commit(&MPI_CNC);
	MPI_Aint lb, extent;
	MPI_Type_get_extent(MPI_Coord, &lb, &extent);
	MPI_Type_create_resized(MPI_CNC, lb, extent, &MPI_CNCR);
	MPI_Type_commit(&MPI_CNCR);
}

int main(int argc, char *argv[]) {

	/* - - - - - Inicialización - - - - - */

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

	ncu = datos.n / npr;
	if((datos.n % npr) > 0) ncu++;
	cuerpos_totales = ncu * npr;

	// Creamos datatype MPI_CNCR

	mpi_datatype_cncr();

	// Se reserva memoria para guardar los datos de los cuerpos

	masas = malloc((sizeof(struct Masas)) * cuerpos_totales);
	if(rank == 0){
		pos = malloc((sizeof(struct Coord)) * cuerpos_totales);
		vel = malloc((sizeof(struct Coord)) * cuerpos_totales);
		acc = malloc((sizeof(struct Coord)) * cuerpos_totales);
	}else{
		vel = malloc((sizeof(struct Coord)) * ncu);
	}

	// Se reserva memoria para los anillos

	p_local = malloc(sizeof(struct Coord) * ncu);
	p_anillo = malloc(sizeof(struct Coord) * ncu);
	a_local = malloc(sizeof(struct Coord) * ncu);
	a_anillo = malloc(sizeof(struct Coord) * ncu);


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

	// Envío/Recepción de las masas (Todos los procesos tienen todas las masas)

	MPI_Bcast(masas, cuerpos_totales, MPI_Masas, 0, MPI_COMM_WORLD);

	// Envío/Recepción de las posiciones y velocidades desde 0 a los demás

	MPI_Scatter(pos, 1, MPI_CNCR, p_local, ncu, MPI_Coord, 0, MPI_COMM_WORLD);
	MPI_Scatter(vel, ncu, MPI_Coord, vel, ncu, MPI_Coord, 0, MPI_COMM_WORLD);


	double inicio, fin; // Tiempo transcurrido

	double t = 0.0;

	int flag = 1;

	MPI_Barrier(MPI_COMM_WORLD);

	inicio = MPI_Wtime();

	// Cálculo de aceleraciones iniciales

	prepararAceleracion();

	MPI_Gather(a_local, ncu, MPI_Coord, acc, ncu, MPI_Coord, 0, MPI_COMM_WORLD);

	for(int pasos = 1; pasos <= datos.tp; pasos++){

		// Recorremos los cuerpos asignados al proceso actual

		for(int i = 0; i < ncu; i++){

			p_local[i].x += vel[i].x * datos.delta;
			p_local[i].y += vel[i].y * datos.delta;
			vel[i].x += a_local[i].x * datos.delta;
			vel[i].y += a_local[i].y * datos.delta;

		}

		// Envío/Recepción hacia/desde todos de las posiciones nuevas

	//	MPI_Allgather(pos, ncu, MPI_Coord, all, ncu, MPI_Coord, MPI_COMM_WORLD);

		// Calculamos aceleraciones

		prepararAceleracion();

		MPI_Gather(p_local, ncu, MPI_Coord, pos, ncu, MPI_Coord, 0, MPI_COMM_WORLD);
		MPI_Gather(vel, ncu, MPI_Coord, vel, ncu, MPI_Coord, 0, MPI_COMM_WORLD);
		MPI_Gather(a_local, ncu, MPI_Coord, acc, ncu, MPI_Coord, 0, MPI_COMM_WORLD);

		// Incrementamos t

		t += datos.delta;

		if(pasos % datos.k == 0){

			// Envío/Recepción de velocidades y aceleraciones

			// MPI_Gather(vel, ncu, MPI_Coord, vel, ncu, MPI_Coord, 0, MPI_COMM_WORLD);
			// MPI_Gather(acc, ncu, MPI_Coord, acc, ncu, MPI_Coord, 0, MPI_COMM_WORLD);

			// Rank 0 imprime cuerpos

			if(rank == 0){
				if(flag){
					 printf("\n            \t%*s\t%*s\t%*s\t%*s\t%*s\t%*s\n", 10, "Posicion(x)", 10, "Posicion(y)", 10, "Velocidad(x)", 10, "Velocidad(y)", 10, "Aceleracion(x)", 10, "Aceleracion(y)");
					 flag = 0;
				 }
				printf("%.2f\n", t);
				imprimirTerminal(0);
			}
		}

	}

	// Calculamos el tiempo de cada proceso

	fin = MPI_Wtime();
	double tiempo_indv = fin - inicio;

	// Obtenemos el tiempo mayor de entre todos

	double tiempo_total;
	MPI_Reduce(&tiempo_indv, &tiempo_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if(rank == 0) printf("\nPrograma ejecutado en %f segundos.\n", tiempo_total);

	/* - - - - Finalización - - - - */

	MPI_Finalize();

	free(masas);
	free(pos);
	free(vel);
	free(acc);
	free(p_local); free(p_anillo);
	free(a_local); free(a_anillo);

	return 0;

}
