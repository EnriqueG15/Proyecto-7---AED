#ifndef MATRIZ_CSR_H
#define MATRIZ_CSR_H

/*
 * Proyecto 7 - Sistema de recomendacion con matrices dispersas
 * Integrante 2: Formato CSR (Compressed Sparse Row - enlaces implicitos)
 * Depende solo de las estructuras Nodo / MatrizDispersa del Integrante 1.
 */

#include <stddef.h>
#include <stdbool.h>
#include "matriz_listas.h"

typedef struct {
    int m;              /* Numero de filas */
    int n;              /* Numero de columnas */
    int nnz;            /* Elementos no nulos */
    float* valores;     /* Arreglo de tamano nnz */
    int* columnas;      /* Arreglo de tamano nnz con las columnas */
    int* ptr_fila;      /* Arreglo de tamano m + 1 con el inicio de cada fila */
} MatrizCSR;

/* ---------------------- Primitivas y conversiones ---------------------- */
MatrizCSR* crearCSR(int m, int n, int nnz);
void destruirCSR(MatrizCSR* A);
MatrizCSR* convertirListasACSR(const MatrizDispersa* M);
float consultarCSR(const MatrizCSR* A, int i, int j);

/* ---------------------- Operaciones matriciales ------------------------ */
MatrizCSR* sumarCSR(const MatrizCSR* A, const MatrizCSR* B);
MatrizCSR* transpuestaCSR(const MatrizCSR* A);
MatrizCSR* multiplicarCSR(const MatrizCSR* A, const MatrizCSR* B);

/* ---------------------- Memoria y utilidades --------------------------- */
size_t memoriaUsadaCSR(const MatrizCSR* A);
void imprimirTripletasCSR(const MatrizCSR* A);
void imprimirEstructuraCSR(const MatrizCSR* A);
void imprimirVistaMatricialCSR(const MatrizCSR* A);
bool sonIgualesCSR(const MatrizCSR* A, const MatrizCSR* B, float tolerancia);

#endif /* MATRIZ_CSR_H */