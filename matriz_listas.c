#include <stdio.h>
#include <stdlib.h>
#include "matriz_listas.h"

static Nodo* obtenerNodo(MatrizDispersa* M, int fila, int col, float calificacion) {
    Nodo* nodo = NULL;
    if (M->avail != NULL) {
        nodo = M->avail;
        M->avail = M->avail->sig_col;
    } else {
        nodo = (Nodo*)malloc(sizeof(Nodo));
        if (!nodo) return NULL;
    }
    nodo->fila = fila;
    nodo->col = col;
    nodo->calificacion = calificacion;
    nodo->sig_col = NULL;
    nodo->sig_fila = NULL;
    return nodo;
}

static void apilarAVAIL(MatrizDispersa* M, Nodo* nodo) {
    nodo->sig_fila = NULL;
    nodo->sig_col = M->avail;
    M->avail = nodo;
}

MatrizDispersa* crearMatriz(int m, int n) {
    if (m <= 0 || n <= 0) return NULL;
    MatrizDispersa* M = (MatrizDispersa*)malloc(sizeof(MatrizDispersa));
    if (!M) return NULL;

    M->m = m;
    M->n = n;
    M->nnz = 0;
    M->avail = NULL;
    M->cabeceras_fila = (Nodo**)calloc((size_t)m, sizeof(Nodo*));
    M->cabeceras_col  = (Nodo**)calloc((size_t)n, sizeof(Nodo*));

    if (!M->cabeceras_fila || !M->cabeceras_col) {
        free(M->cabeceras_fila);
        free(M->cabeceras_col);
        free(M);
        return NULL;
    }
    return M;
}

float consultarCalificacion(const MatrizDispersa* M, int i, int j) {
    if (!M || i < 0 || i >= M->m || j < 0 || j >= M->n) return 0.0f;
    const Nodo* curr = M->cabeceras_fila[i];
    while (curr != NULL && curr->col < j) curr = curr->sig_col;
    if (curr != NULL && curr->col == j) return curr->calificacion;
    return 0.0f;
}

void eliminarCalificacion(MatrizDispersa* M, int i, int j) {
    if (!M || i < 0 || i >= M->m || j < 0 || j >= M->n) return;

    Nodo* prev_f = NULL;
    Nodo* curr_f = M->cabeceras_fila[i];
    while (curr_f != NULL && curr_f->col < j) {
        prev_f = curr_f;
        curr_f = curr_f->sig_col;
    }
    if (curr_f == NULL || curr_f->col != j) return;

    Nodo* prev_c = NULL;
    Nodo* curr_c = M->cabeceras_col[j];
    while (curr_c != NULL && curr_c->fila < i) {
        prev_c = curr_c;
        curr_c = curr_c->sig_fila;
    }

    if (prev_f == NULL) M->cabeceras_fila[i] = curr_f->sig_col;
    else prev_f->sig_col = curr_f->sig_col;

    if (prev_c == NULL) M->cabeceras_col[j] = curr_c->sig_fila;
    else prev_c->sig_fila = curr_c->sig_fila;

    apilarAVAIL(M, curr_f);
    M->nnz--;
}

void insertarCalificacion(MatrizDispersa* M, int i, int j, float val) {
    if (!M || i < 0 || i >= M->m || j < 0 || j >= M->n) return;
    if (val == 0.0f) {
        eliminarCalificacion(M, i, j);
        return;
    }

    Nodo* prev_f = NULL;
    Nodo* curr_f = M->cabeceras_fila[i];
    while (curr_f != NULL && curr_f->col < j) {
        prev_f = curr_f;
        curr_f = curr_f->sig_col;
    }
    if (curr_f != NULL && curr_f->col == j) {
        curr_f->calificacion = val;
        return;
    }

    Nodo* prev_c = NULL;
    Nodo* curr_c = M->cabeceras_col[j];
    while (curr_c != NULL && curr_c->fila < i) {
        prev_c = curr_c;
        curr_c = curr_c->sig_fila;
    }

    Nodo* nuevo = obtenerNodo(M, i, j, val);
    if (!nuevo) return;

    nuevo->sig_col = curr_f;
    if (prev_f == NULL) M->cabeceras_fila[i] = nuevo;
    else prev_f->sig_col = nuevo;

    nuevo->sig_fila = curr_c;
    if (prev_c == NULL) M->cabeceras_col[j] = nuevo;
    else prev_c->sig_fila = nuevo;

    M->nnz++;
}

void destruirMatrizDispersa(MatrizDispersa* M) {
    if (!M) return;
    for (int i = 0; i < M->m; i++) {
        Nodo* curr = M->cabeceras_fila[i];
        while (curr != NULL) {
            Nodo* sig = curr->sig_col;
            free(curr);
            curr = sig;
        }
    }
    Nodo* av = M->avail;
    while (av != NULL) {
        Nodo* sig = av->sig_col;
        free(av);
        av = sig;
    }
    free(M->cabeceras_fila);
    free(M->cabeceras_col);
    free(M);
}

void imprimirReporteMemoria(const MatrizDispersa* M) {
    if (!M) return;
    int n_avail = 0;
    for (const Nodo* p = M->avail; p != NULL; p = p->sig_col) n_avail++;

    size_t b_nodos = (size_t)M->nnz * sizeof(Nodo);
    size_t b_cab   = ((size_t)M->m + (size_t)M->n) * sizeof(Nodo*);
    size_t b_avail = (size_t)n_avail * sizeof(Nodo);
    size_t b_disp  = sizeof(MatrizDispersa) + b_nodos + b_cab + b_avail;
    size_t b_densa = (size_t)M->m * (size_t)M->n * sizeof(float) + (size_t)M->m * sizeof(float*);
    double celdas  = (double)M->m * (double)M->n;

    printf("\n=== REPORTE DE MEMORIA (LISTAS CRUZADAS VS DENSA) ===\n");
    printf("Dimensiones: %d x %d | nnz: %d | Dispersion: %.4f%%\n",
           M->m, M->n, M->nnz, (1.0 - (double)M->nnz / celdas) * 100.0);
    printf("Listas Cruzadas: %zu bytes (%.2f MB) [nodos: %zu, cabeceras: %zu, AVAIL(%d): %zu]\n",
           b_disp, (double)b_disp / (1024.0 * 1024.0), b_nodos, b_cab, n_avail, b_avail);
    printf("Matriz Densa   : %zu bytes (%.2f MB)\n",
           b_densa, (double)b_densa / (1024.0 * 1024.0));
    if (b_densa > 0) {
        printf("Ahorro neto    : %.2f%%\n", (1.0 - (double)b_disp / (double)b_densa) * 100.0);
    }
}