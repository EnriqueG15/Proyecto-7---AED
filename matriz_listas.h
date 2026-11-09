#ifndef MATRIZ_LISTAS_H
#define MATRIZ_LISTAS_H

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

#endif