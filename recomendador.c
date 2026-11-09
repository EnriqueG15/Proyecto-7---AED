/*
 * recomendador.c  -  Integrante 3  (lenguaje C)
 * Todas las estructuras con typedef struct + malloc/calloc/free.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "recomendador.h"

/* ====================== Orden de ranking ======================
 * "Mejor" = mayor puntaje; si empatan, menor id (resultados deterministas).
 * HEAP_MAX: la raiz es el mejor.  HEAP_MIN: la raiz es el PEOR elemento. */
static bool rankMejor(int idA, float pA, int idB, float pB) {
    if (pA != pB) return pA > pB;
    return idA < idB;
}

static bool antes(const TColaPrioridad* C, const ParHeap* a, const ParHeap* b) {
    if (C->modo == HEAP_MAX) return rankMejor(a->id, a->puntaje, b->id, b->puntaje);
    return rankMejor(b->id, b->puntaje, a->id, a->puntaje);
}

static void intercambiar(ParHeap* a, ParHeap* b) {
    ParHeap t = *a; *a = *b; *b = t;
}

/* flotar: sube el elemento i hacia la raiz mientras "antes" que su padre */
static void flotar(TColaPrioridad* C, int i) {
    while (i > 0) {
        int p = (i - 1) / 2;
        if (antes(C, &C->arreglo[i], &C->arreglo[p])) {
            intercambiar(&C->arreglo[i], &C->arreglo[p]);
            i = p;
        } else break;
    }
}

/* hundir: baja el elemento i intercambiando con el hijo "optimo" */
static void hundir(TColaPrioridad* C, int i) {
    for (;;) {
        int l = 2 * i + 1, r = l + 1, o = i;
        if (l < C->tamano && antes(C, &C->arreglo[l], &C->arreglo[o])) o = l;
        if (r < C->tamano && antes(C, &C->arreglo[r], &C->arreglo[o])) o = r;
        if (o == i) break;
        intercambiar(&C->arreglo[i], &C->arreglo[o]);
        i = o;
    }
}

/* ====================== TAD ColaPrioridad ====================== */

TColaPrioridad* crearColaPrioridad(int capacidad, ModoHeap modo) {
    if (capacidad <= 0) return NULL;
    TColaPrioridad* C = (TColaPrioridad*)malloc(sizeof(TColaPrioridad));
    if (!C) return NULL;
    C->arreglo = (ParHeap*)malloc((size_t)capacidad * sizeof(ParHeap));
    if (!C->arreglo) { free(C); return NULL; }
    C->tamano = 0;
    C->capacidad = capacidad;
    C->modo = modo;
    return C;
}

void destruirColaPrioridad(TColaPrioridad* C) {
    if (!C) return;
    free(C->arreglo);
    free(C);
}

bool esVacia(const TColaPrioridad* C) { return C->tamano == 0; }

bool insertarHeap(TColaPrioridad* C, int id, float puntaje) {
    if (!C || C->tamano >= C->capacidad) return false;
    C->arreglo[C->tamano].id = id;
    C->arreglo[C->tamano].puntaje = puntaje;
    C->tamano++;
    flotar(C, C->tamano - 1);
    return true;
}

ParHeap extraerOptimo(TColaPrioridad* C) {
    ParHeap vacio = { -1, 0.0f };
    if (!C || C->tamano == 0) return vacio;
    ParHeap raiz = C->arreglo[0];
    C->tamano--;
    if (C->tamano > 0) {
        C->arreglo[0] = C->arreglo[C->tamano];
        hundir(C, 0);
    }
    return raiz;
}

ParHeap verOptimo(const TColaPrioridad* C) {
    ParHeap vacio = { -1, 0.0f };
    if (!C || C->tamano == 0) return vacio;
    return C->arreglo[0];
}

void reemplazarOptimo(TColaPrioridad* C, int id, float puntaje) {
    if (!C || C->tamano == 0) return;
    C->arreglo[0].id = id;
    C->arreglo[0].puntaje = puntaje;
    hundir(C, 0);
}

/* ====================== Dataset ====================== */

int cargarMovieLens(MatrizDispersa* M, const char* ruta) {
    FILE* f = fopen(ruta, "r");
    if (!f) return -1;

    char linea[256];
    int leidas = 0;
    while (fgets(linea, sizeof(linea), f)) {
        int u, i;
        float r;
        /* acepta tabuladores, espacios o comas; ignora cabeceras */
        if (sscanf(linea, "%d%*[,\t ;]%d%*[,\t ;]%f", &u, &i, &r) != 3) continue;
        u -= 1;  /* los ids de MovieLens empiezan en 1 */
        i -= 1;
        if (u < 0 || u >= M->m || i < 0 || i >= M->n || r <= 0.0f) continue;
        insertarCalificacion(M, u, i, r);
        leidas++;
    }
    fclose(f);
    return leidas;
}

static unsigned int siguienteAleatorio(unsigned int* s) {  /* xorshift32, portable */
    unsigned int x = *s;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *s = x;
    return x;
}

int generarSintetico(MatrizDispersa* M, int nnzObjetivo, unsigned int semilla) {
    if (!M || nnzObjetivo <= 0) return M ? M->nnz : 0;
    double celdas = (double)M->m * (double)M->n;
    long long maxPermitido = (long long)(0.05 * celdas);
    long long objetivo = nnzObjetivo;
    if (objetivo > maxPermitido) {
        fprintf(stderr, "Aviso: se limita nnz a %lld para mantener >95%% de celdas vacias.\n",
                maxPermitido);
        objetivo = maxPermitido;
    }
    if (semilla == 0) semilla = 2463534242u;

    long long intentos = 0, limite = objetivo * 50 + 1000;
    while (M->nnz < objetivo && intentos < limite) {
        int u = (int)(siguienteAleatorio(&semilla) % (unsigned int)M->m);
        int i = (int)(siguienteAleatorio(&semilla) % (unsigned int)M->n);
        float v = (float)(1 + siguienteAleatorio(&semilla) % 5u);
        insertarCalificacion(M, u, i, v);  /* si ya existia, solo actualiza */
        intentos++;
    }
    return M->nnz;
}

/* ====================== Similitud coseno ====================== */

float* calcularNormas(const MatrizCSR* A) {
    float* normas = (float*)calloc((size_t)A->m, sizeof(float));
    if (!normas) return NULL;
    for (int u = 0; u < A->m; u++) {
        double s = 0.0;
        for (int p = A->ptr_fila[u]; p < A->ptr_fila[u + 1]; p++)
            s += (double)A->valores[p] * (double)A->valores[p];
        normas[u] = (float)sqrt(s);
    }
    return normas;
}

float similitudCoseno(const MatrizCSR* A, const float* normas, int u, int v) {
    if (!A || !normas || u < 0 || u >= A->m || v < 0 || v >= A->m) return 0.0f;
    double den = (double)normas[u] * (double)normas[v];
    if (den == 0.0) return 0.0f;

    int i = A->ptr_fila[u], fi = A->ptr_fila[u + 1];
    int j = A->ptr_fila[v], fj = A->ptr_fila[v + 1];
    double dot = 0.0;
    while (i < fi && j < fj) {           /* interseccion de secuencias ordenadas */
        int cu = A->columnas[i], cv = A->columnas[j];
        if (cu == cv) { dot += (double)A->valores[i] * (double)A->valores[j]; i++; j++; }
        else if (cu < cv) i++;
        else j++;
    }
    return (float)(dot / den);
}

/* ====================== Vecinos y Top-k ====================== */

int usuariosSimilares(const MatrizCSR* A, const float* normas, int u, int nVecinos, ParHeap* salida) {
    if (!A || !normas || u < 0 || u >= A->m || nVecinos <= 0 || normas[u] == 0.0f) return 0;

    /* min-heap de tamano nVecinos: la raiz es el peor de los mejores */
    TColaPrioridad* H = crearColaPrioridad(nVecinos, HEAP_MIN);
    if (!H) return 0;

    for (int v = 0; v < A->m; v++) {
        if (v == u) continue;
        float s = similitudCoseno(A, normas, u, v);
        if (s <= 0.0f) continue;
        if (H->tamano < nVecinos) {
            insertarHeap(H, v, s);
        } else {
            ParHeap peor = verOptimo(H);
            if (rankMejor(v, s, peor.id, peor.puntaje)) reemplazarOptimo(H, v, s);
        }
    }
    int cnt = H->tamano;
    for (int p = cnt - 1; p >= 0; p--) salida[p] = extraerOptimo(H);  /* peor -> mejor, llenando desde atras */
    destruirColaPrioridad(H);
    return cnt;
}

int recomendarTopK(const MatrizCSR* A, const float* normas, int u, int k,
                   int nVecinos, Recomendacion* salida) {
    if (!A || !normas || u < 0 || u >= A->m || k <= 0 || nVecinos <= 0) return 0;

    ParHeap* vec = (ParHeap*)malloc((size_t)nVecinos * sizeof(ParHeap));
    if (!vec) return 0;
    int nv = usuariosSimilares(A, normas, u, nVecinos, vec);
    if (nv == 0) { free(vec); return 0; }

    int n = A->n;
    double* num = (double*)calloc((size_t)n, sizeof(double));
    double* den = (double*)calloc((size_t)n, sizeof(double));
    char*   visto = (char*)calloc((size_t)n, sizeof(char));
    int*    cand = (int*)malloc((size_t)n * sizeof(int));
    TColaPrioridad* H = crearColaPrioridad(k, HEAP_MIN);
    if (!num || !den || !visto || !cand || !H) {
        free(vec); free(num); free(den); free(visto); free(cand); destruirColaPrioridad(H);
        return 0;
    }

    /* items ya consumidos por u: no se recomiendan */
    for (int p = A->ptr_fila[u]; p < A->ptr_fila[u + 1]; p++) visto[A->columnas[p]] = 1;

    /* acumular sum(sim*rating) y sum(sim) por item candidato */
    int nc = 0;
    for (int j = 0; j < nv; j++) {
        int v = vec[j].id;
        double s = vec[j].puntaje;
        for (int p = A->ptr_fila[v]; p < A->ptr_fila[v + 1]; p++) {
            int it = A->columnas[p];
            if (visto[it]) continue;
            if (den[it] == 0.0) cand[nc++] = it;
            num[it] += s * (double)A->valores[p];
            den[it] += s;
        }
    }

    /* Top-k con min-heap acotado a k: O(I log k) */
    for (int c = 0; c < nc; c++) {
        int it = cand[c];
        float score = (float)(num[it] / den[it]);
        if (H->tamano < k) {
            insertarHeap(H, it, score);
        } else {
            ParHeap peor = verOptimo(H);
            if (rankMejor(it, score, peor.id, peor.puntaje)) reemplazarOptimo(H, it, score);
        }
    }

    int cnt = H->tamano;
    for (int p = cnt - 1; p >= 0; p--) {
        ParHeap e = extraerOptimo(H);
        salida[p].id_item = e.id;
        salida[p].puntaje = e.puntaje;
    }

    destruirColaPrioridad(H);
    free(vec); free(num); free(den); free(visto); free(cand);
    return cnt;
}