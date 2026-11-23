#ifndef MATRIZ_LISTAS_H
#define MATRIZ_LISTAS_H

#include <stdbool.h>

typedef struct Nodo {
    int fila;
    int col;
    float calificacion;
    struct Nodo* sig_col;
    struct Nodo* sig_fila;
} Nodo;

typedef struct MatrizDispersa {
    int m;
    int n;
    int nnz;
    Nodo** cabeceras_fila;
    Nodo** cabeceras_col;
    Nodo* avail;
} MatrizDispersa;

MatrizDispersa* crearMatriz(int m, int n);
void destruirMatrizDispersa(MatrizDispersa* M);
void insertarCalificacion(MatrizDispersa* M, int i, int j, float val);
float consultarCalificacion(const MatrizDispersa* M, int i, int j);
void eliminarCalificacion(MatrizDispersa* M, int i, int j);
void imprimirReporteMemoria(const MatrizDispersa* M);

#endif