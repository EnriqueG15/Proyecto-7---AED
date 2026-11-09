#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "matriz_csr.h"

#define MAX_IMPRESION 100

static float absf(float x) { return x < 0.0f ? -x : x; }

static int cmpInt(const void* a, const void* b) {
    int x = *(const int*)a, y = *(const int*)b;
    return (x > y) - (x < y);
}

/* ============================ PRIMITIVAS ============================ */

MatrizCSR* crearCSR(int m, int n, int nnz) {
    if (m <= 0 || n <= 0 || nnz < 0) return NULL;

    MatrizCSR* A = (MatrizCSR*)malloc(sizeof(MatrizCSR));
    if (A == NULL) return NULL;

    size_t cap = (nnz > 0) ? (size_t)nnz : 1;   /* evita malloc(0) */
    A->m = m;
    A->n = n;
    A->nnz = nnz;
    A->valores  = (float*)malloc(cap * sizeof(float));
    A->columnas = (int*)malloc(cap * sizeof(int));
    A->ptr_fila = (int*)calloc((size_t)m + 1, sizeof(int));

    if (A->valores == NULL || A->columnas == NULL || A->ptr_fila == NULL) {
        free(A->valores);
        free(A->columnas);
        free(A->ptr_fila);
        free(A);
        return NULL;
    }
    return A;
}

void destruirCSR(MatrizCSR* A) {
    if (A == NULL) return;
    free(A->valores);
    free(A->columnas);
    free(A->ptr_fila);
    free(A);
}

/* Recorta valores/columnas al tamano exacto nnz */
static void ajustarCapacidadCSR(MatrizCSR* A) {
    size_t cap = (A->nnz > 0) ? (size_t)A->nnz : 1;
    float* v = (float*)realloc(A->valores, cap * sizeof(float));
    int* c = (int*)realloc(A->columnas, cap * sizeof(int));
    if (v != NULL) A->valores = v;
    if (c != NULL) A->columnas = c;
}

/* Aumenta la capacidad (duplicando) hasta poder alojar 'requerido' elementos */
static bool crecerCSR(MatrizCSR* C, int* cap, int requerido) {
    long long nueva = (*cap > 0) ? *cap : 1;
    while (nueva < requerido) nueva *= 2;
    if (nueva > 2147483647LL) nueva = requerido;

    float* v = (float*)realloc(C->valores, (size_t)nueva * sizeof(float));
    if (v == NULL) return false;
    C->valores = v;

    int* c = (int*)realloc(C->columnas, (size_t)nueva * sizeof(int));
    if (c == NULL) return false;
    C->columnas = c;

    *cap = (int)nueva;
    return true;
}

/* ============================ CONVERSION ============================ */
/*
 * Recorre cada lista horizontal (cabeceras_fila[i] -> sig_col).
 * Pasada 1: cuenta los nodos reales (no depende de que M->nnz este al dia).
 * Pasada 2: llena columnas/valores y registra ptr_fila.
 * Complejidad total O(m + nnz).
 */
MatrizCSR* convertirListasACSR(const MatrizDispersa* M) {
    if (M == NULL) return NULL;

    int cuenta = 0;
    for (int i = 0; i < M->m; i++)
        for (const Nodo* p = M->cabeceras_fila[i]; p != NULL; p = p->sig_col)
            cuenta++;

    MatrizCSR* A = crearCSR(M->m, M->n, cuenta);
    if (A == NULL) return NULL;

    int k = 0;
    for (int i = 0; i < M->m; i++) {
        A->ptr_fila[i] = k;
        for (const Nodo* p = M->cabeceras_fila[i]; p != NULL; p = p->sig_col) {
            A->columnas[k] = p->col;   /* las listas del Int. 1 estan ordenadas por columna */
            A->valores[k]  = p->calificacion;
            k++;
        }
    }
    A->ptr_fila[M->m] = k;
    return A;
}

/* ============================ CONSULTA ============================ */
/* Busqueda binaria en el segmento [ptr_fila[i], ptr_fila[i+1]-1]: O(log nnz_fila) */
float consultarCSR(const MatrizCSR* A, int i, int j) {
    if (A == NULL || i < 0 || i >= A->m || j < 0 || j >= A->n) return 0.0f;

    int lo = A->ptr_fila[i];
    int hi = A->ptr_fila[i + 1] - 1;
    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;
        int c = A->columnas[mid];
        if (c == j) return A->valores[mid];
        if (c < j) lo = mid + 1;
        else       hi = mid - 1;
    }
    return 0.0f;
}

/* ============================ SUMA ============================ */
/* Merge de 2 punteros por fila: O(m + nnz_A + nnz_B). Omite resultados 0.0f */
MatrizCSR* sumarCSR(const MatrizCSR* A, const MatrizCSR* B) {
    if (A == NULL || B == NULL || A->m != B->m || A->n != B->n) return NULL;

    MatrizCSR* C = crearCSR(A->m, A->n, A->nnz + B->nnz);  /* cota superior */
    if (C == NULL) return NULL;

    int k = 0;
    for (int i = 0; i < A->m; i++) {
        C->ptr_fila[i] = k;
        int pa = A->ptr_fila[i], fa = A->ptr_fila[i + 1];
        int pb = B->ptr_fila[i], fb = B->ptr_fila[i + 1];

        while (pa < fa && pb < fb) {
            int ca = A->columnas[pa];
            int cb = B->columnas[pb];
            if (ca < cb) {
                C->columnas[k] = ca; C->valores[k] = A->valores[pa]; k++; pa++;
            } else if (cb < ca) {
                C->columnas[k] = cb; C->valores[k] = B->valores[pb]; k++; pb++;
            } else {
                float s = A->valores[pa] + B->valores[pb];
                if (s != 0.0f) { C->columnas[k] = ca; C->valores[k] = s; k++; }
                pa++; pb++;
            }
        }
        while (pa < fa) { C->columnas[k] = A->columnas[pa]; C->valores[k] = A->valores[pa]; k++; pa++; }
        while (pb < fb) { C->columnas[k] = B->columnas[pb]; C->valores[k] = B->valores[pb]; k++; pb++; }
    }
    C->ptr_fila[A->m] = k;
    C->nnz = k;
    ajustarCapacidadCSR(C);
    return C;
}

/* ============================ TRANSPUESTA ============================ */
/* 3 pasadas tipo Counting Sort: O(m + n + nnz) */
MatrizCSR* transpuestaCSR(const MatrizCSR* A) {
    if (A == NULL) return NULL;

    MatrizCSR* T = crearCSR(A->n, A->m, A->nnz);   /* ptr_fila (n+1) inicia en 0 */
    if (T == NULL) return NULL;

    /* (a) Conteo de frecuencias: la columna j de A sera la fila j de A^T.
           Se cuenta en j+1 para que la suma prefija deje el inicio en j. */
    for (int k = 0; k < A->nnz; k++)
        T->ptr_fila[A->columnas[k] + 1]++;

    /* (b) Sumas prefijas acumuladas */
    for (int j = 0; j < A->n; j++)
        T->ptr_fila[j + 1] += T->ptr_fila[j];

    /* (c) Distribucion directa con cursores */
    int* cursor = (int*)malloc((size_t)A->n * sizeof(int));
    if (cursor == NULL) { destruirCSR(T); return NULL; }
    for (int j = 0; j < A->n; j++) cursor[j] = T->ptr_fila[j];

    for (int i = 0; i < A->m; i++) {
        for (int k = A->ptr_fila[i]; k < A->ptr_fila[i + 1]; k++) {
            int j = A->columnas[k];
            int pos = cursor[j]++;
            T->columnas[pos] = i;           /* i crece => columnas de A^T ordenadas */
            T->valores[pos]  = A->valores[k];
        }
    }
    free(cursor);
    return T;
}

/* ============================ PRODUCTO ============================ */
/*
 * Algoritmo de Gustavson (acumulador disperso por fila):
 *   fila i de C = suma_k  A[i][k] * (fila k de B)
 * - acc[n]      : acumulador de valores.
 * - marcador[n] : marcador[j] == i  <=> la columna j ya aparecio en la fila i
 *                 (evita reinicializar acc completo en cada fila).
 * - cols_fila   : lista de columnas tocadas en la fila actual.
 * Para dejar las columnas ordenadas se elige lo mas barato:
 *   ordenar cols_fila (cnt*log cnt) o barrer marcador en O(n).
 */
MatrizCSR* multiplicarCSR(const MatrizCSR* A, const MatrizCSR* B) {
    if (A == NULL || B == NULL || A->n != B->m) return NULL;

    int m = A->m, n = B->n;
    int cap = A->nnz + B->nnz;
    if (cap < 16) cap = 16;

    MatrizCSR* C     = crearCSR(m, n, cap);
    float* acc       = (float*)calloc((size_t)n, sizeof(float));
    int* marcador    = (int*)malloc((size_t)n * sizeof(int));
    int* cols_fila   = (int*)malloc((size_t)n * sizeof(int));

    if (C == NULL || acc == NULL || marcador == NULL || cols_fila == NULL) {
        destruirCSR(C); free(acc); free(marcador); free(cols_fila);
        return NULL;
    }
    for (int j = 0; j < n; j++) marcador[j] = -1;

    int k = 0;
    for (int i = 0; i < m; i++) {
        C->ptr_fila[i] = k;
        int cnt = 0;

        for (int pa = A->ptr_fila[i]; pa < A->ptr_fila[i + 1]; pa++) {
            int kk  = A->columnas[pa];
            float a = A->valores[pa];
            for (int pb = B->ptr_fila[kk]; pb < B->ptr_fila[kk + 1]; pb++) {
                int j = B->columnas[pb];
                if (marcador[j] != i) {
                    marcador[j] = i;
                    cols_fila[cnt++] = j;
                    acc[j] = a * B->valores[pb];
                } else {
                    acc[j] += a * B->valores[pb];
                }
            }
        }
        if (cnt == 0) continue;

        /* Orden ascendente de columnas */
        int lg = 0;
        while ((1 << lg) < cnt && lg < 30) lg++;
        if ((long long)cnt * (lg + 1) < (long long)n) {
            qsort(cols_fila, (size_t)cnt, sizeof(int), cmpInt);
        } else {
            int c2 = 0;
            for (int j = 0; j < n; j++)
                if (marcador[j] == i) cols_fila[c2++] = j;
        }

        if (k + cnt > cap && !crecerCSR(C, &cap, k + cnt)) {
            destruirCSR(C); free(acc); free(marcador); free(cols_fila);
            return NULL;
        }
        for (int t = 0; t < cnt; t++) {
            int j = cols_fila[t];
            if (acc[j] != 0.0f) {            /* descarta cancelaciones exactas */
                C->columnas[k] = j;
                C->valores[k]  = acc[j];
                k++;
            }
        }
    }
    C->ptr_fila[m] = k;
    C->nnz = k;
    ajustarCapacidadCSR(C);

    free(acc);
    free(marcador);
    free(cols_fila);
    return C;
}

/* ============================ MEMORIA ============================ */
size_t memoriaUsadaCSR(const MatrizCSR* A) {
    if (A == NULL) return 0;
    return sizeof(MatrizCSR)
         + (size_t)A->nnz * sizeof(float)
         + (size_t)A->nnz * sizeof(int)
         + (size_t)(A->m + 1) * sizeof(int);
}

/* ============================ IMPRESION ============================ */
void imprimirTripletasCSR(const MatrizCSR* A) {
    if (A == NULL) { printf("  (matriz no disponible)\n"); return; }
    if (A->nnz == 0) { printf("  (matriz nula: sin elementos no nulos)\n"); return; }

    int c = 0;
    for (int i = 0; i < A->m && c < MAX_IMPRESION; i++)
        for (int k = A->ptr_fila[i]; k < A->ptr_fila[i + 1] && c < MAX_IMPRESION; k++, c++)
            printf("  (%d, %d, %.2f)\n", i, A->columnas[k], A->valores[k]);
    if (A->nnz > c) printf("  ... (%d tripletas mas)\n", A->nnz - c);
}

static void imprimirArregloInt(const char* nombre, const int* a, int len) {
    int lim = len < MAX_IMPRESION ? len : MAX_IMPRESION;
    printf("  %-9s[%d] = [", nombre, len);
    for (int t = 0; t < lim; t++) printf("%s%d", t ? ", " : "", a[t]);
    if (len > lim) printf(", ... (+%d)", len - lim);
    printf("]\n");
}

static void imprimirArregloFloat(const char* nombre, const float* a, int len) {
    int lim = len < MAX_IMPRESION ? len : MAX_IMPRESION;
    printf("  %-9s[%d] = [", nombre, len);
    for (int t = 0; t < lim; t++) printf("%s%.2f", t ? ", " : "", a[t]);
    if (len > lim) printf(", ... (+%d)", len - lim);
    printf("]\n");
}

void imprimirEstructuraCSR(const MatrizCSR* A) {
    if (A == NULL) { printf("  (matriz no disponible)\n"); return; }
    imprimirArregloInt("ptr_fila", A->ptr_fila, A->m + 1);
    imprimirArregloInt("columnas", A->columnas, A->nnz);
    imprimirArregloFloat("valores", A->valores, A->nnz);
}

void imprimirVistaMatricialCSR(const MatrizCSR* A) {
    if (A == NULL) return;
    if (A->m > 15 || A->n > 15) {
        printf("  (vista matricial omitida: %d x %d supera 15 x 15)\n", A->m, A->n);
        return;
    }
    printf("       ");
    for (int j = 0; j < A->n; j++) printf("%8d", j);
    printf("\n");
    for (int i = 0; i < A->m; i++) {
        printf("  %3d |", i);
        int k = A->ptr_fila[i];
        for (int j = 0; j < A->n; j++) {
            if (k < A->ptr_fila[i + 1] && A->columnas[k] == j) {
                printf("%8.2f", A->valores[k]);
                k++;
            } else {
                printf("%8.2f", 0.0f);
            }
        }
        printf("\n");
    }
}

bool sonIgualesCSR(const MatrizCSR* A, const MatrizCSR* B, float tolerancia) {
    if (A == NULL || B == NULL) return false;
    if (A->m != B->m || A->n != B->n || A->nnz != B->nnz) return false;
    for (int i = 0; i <= A->m; i++)
        if (A->ptr_fila[i] != B->ptr_fila[i]) return false;
    for (int k = 0; k < A->nnz; k++) {
        if (A->columnas[k] != B->columnas[k]) return false;
        if (absf(A->valores[k] - B->valores[k]) > tolerancia) return false;
    }
    return true;
}