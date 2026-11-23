#ifndef RECOMENDADOR_H
#define RECOMENDADOR_H

#include <stddef.h>
#include <stdbool.h>
#include "matriz_listas.h"
#include "matriz_csr.h"

typedef enum { HEAP_MAX = 0, HEAP_MIN = 1 } ModoHeap;

typedef struct {
    int   id;
    float puntaje;
} ParHeap;

typedef struct {
    ParHeap* arreglo;
    int      tamano;
    int      capacidad;
    ModoHeap modo;
} TColaPrioridad;

TColaPrioridad* crearColaPrioridad(int capacidad, ModoHeap modo);
void            destruirColaPrioridad(TColaPrioridad* C);
bool    esVacia(const TColaPrioridad* C);
bool    insertarHeap(TColaPrioridad* C, int id, float puntaje);
ParHeap extraerOptimo(TColaPrioridad* C);
ParHeap verOptimo(const TColaPrioridad* C);
void    reemplazarOptimo(TColaPrioridad* C, int id, float puntaje);

int cargarMovieLens(MatrizDispersa* M, const char* ruta);
int generarSintetico(MatrizDispersa* M, int nnzObjetivo, unsigned int semilla);

typedef struct {
    int   id_item;
    float puntaje;
} Recomendacion;

float* calcularNormas(const MatrizCSR* A);
float similitudCoseno(const MatrizCSR* A, const float* normas, int u, int v);
int usuariosSimilares(const MatrizCSR* A, const float* normas, int u, int nVecinos, ParHeap* salida);
int recomendarTopK(const MatrizCSR* A, const float* normas, int u, int k,
                   int nVecinos, Recomendacion* salida);

#endif