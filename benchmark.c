#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "matriz_listas.h"
#include "matriz_csr.h"
#include "benchmark.h"

static unsigned int rand32(void) {
    return ((unsigned int)rand() << 15) ^ (unsigned int)rand();
}

static MatrizDispersa* crearListasBench(int m, int n) {
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

static bool insertarListasBench(MatrizDispersa* M, int i, int j, float v) {
    Nodo* ant = NULL;
    Nodo* p = M->cabeceras_fila[i];
    while (p && p->col < j) { ant = p; p = p->sig_col; }
    if (p && p->col == j) return false;

    Nodo* q = (Nodo*)malloc(sizeof(Nodo));
    if (!q) return false;
    q->fila = i; q->col = j; q->calificacion = v;
    q->sig_col = p;
    if (ant) ant->sig_col = q; else M->cabeceras_fila[i] = q;

    Nodo* antc = NULL;
    Nodo* r = M->cabeceras_col[j];
    while (r && r->fila < i) { antc = r; r = r->sig_fila; }
    q->sig_fila = r;
    if (antc) antc->sig_fila = q; else M->cabeceras_col[j] = q;

    M->nnz++;
    return true;
}

static float consultarListasBench(const MatrizDispersa* M, int i, int j) {
    const Nodo* p = M->cabeceras_fila[i];
    while (p && p->col < j) p = p->sig_col;
    return (p && p->col == j) ? p->calificacion : 0.0f;
}

static void liberarListasBench(MatrizDispersa* M) {
    if (!M) return;
    for (int i = 0; i < M->m; i++) {
        Nodo* p = M->cabeceras_fila[i];
        while (p) { Nodo* s = p->sig_col; free(p); p = s; }
    }
    free(M->cabeceras_fila);
    free(M->cabeceras_col);
    free(M);
}

static double msTranscurridos(clock_t ini, clock_t fin) {
    return (double)(fin - ini) * 1000.0 / CLOCKS_PER_SEC;
}

void ejecutarBenchmarkPersonalizado(int m, int n, int nnz, unsigned int semilla, int consultas) {
    long long max_celdas = (long long)m * n;
    if (m <= 0 || n <= 0 || nnz < 0 || nnz > max_celdas) {
        printf("Parametros invalidos.\n");
        return;
    }

    srand(semilla);
    MatrizDispersa* M = crearListasBench(m, n);
    if (!M) { printf("Error de memoria.\n"); return; }

    while (M->nnz < nnz) {
        int i = (int)(rand32() % (unsigned int)m);
        int j = (int)(rand32() % (unsigned int)n);
        float v = (float)((rand() % 10) + 1) * 0.5f;
        insertarListasBench(M, i, j, v);
    }

    int* qf = (int*)malloc((size_t)consultas * sizeof(int));
    int* qc = (int*)malloc((size_t)consultas * sizeof(int));
    for (int k = 0; k < consultas; k++) {
        qf[k] = (int)(rand32() % (unsigned int)m);
        qc[k] = (int)(rand32() % (unsigned int)n);
    }

    clock_t t0 = clock();
    MatrizCSR* A = convertirListasACSR(M);
    clock_t t1 = clock();
    double t_conv = msTranscurridos(t0, t1);

    volatile float sink = 0.0f;
    t0 = clock();
    for (int k = 0; k < consultas; k++) sink += consultarListasBench(M, qf[k], qc[k]);
    t1 = clock();
    double t_cons_listas = msTranscurridos(t0, t1);

    t0 = clock();
    for (int k = 0; k < consultas; k++) sink += consultarCSR(A, qf[k], qc[k]);
    t1 = clock();
    double t_cons_csr = msTranscurridos(t0, t1);
    (void)sink;

    t0 = clock();
    MatrizCSR* S = sumarCSR(A, A);
    t1 = clock();
    double t_suma = msTranscurridos(t0, t1);

    t0 = clock();
    MatrizCSR* T = transpuestaCSR(A);
    t1 = clock();
    double t_trans = msTranscurridos(t0, t1);

    t0 = clock();
    MatrizCSR* P = multiplicarCSR(A, T);
    t1 = clock();
    double t_prod = msTranscurridos(t0, t1);

    size_t mem_densa  = (size_t)m * (size_t)n * sizeof(float) + (size_t)m * sizeof(float*);
    size_t mem_listas = sizeof(MatrizDispersa) + (size_t)(m + n) * sizeof(Nodo*) + (size_t)M->nnz * sizeof(Nodo);
    size_t mem_csr    = memoriaUsadaCSR(A);
    double disp = (1.0 - (double)nnz / (double)max_celdas) * 100.0;

    printf("\n=== BENCHMARK: %d x %d | nnz = %d (Dispersion: %.2f%%) | Semilla: %u ===\n",
           m, n, nnz, disp, semilla);
    printf("MEMORIA:\n");
    printf("  Densa        : %12zu bytes (%.2f MB)\n", mem_densa,  mem_densa  / (1024.0 * 1024.0));
    printf("  Listas Cruz. : %12zu bytes (%.2f MB) | Ahorro vs Densa: %.2f%%\n",
           mem_listas, mem_listas / (1024.0 * 1024.0), (1.0 - (double)mem_listas / mem_densa) * 100.0);
    printf("  Formato CSR  : %12zu bytes (%.2f MB) | Ahorro vs Densa: %.2f%%\n",
           mem_csr,    mem_csr    / (1024.0 * 1024.0), (1.0 - (double)mem_csr / mem_densa) * 100.0);

    printf("TIEMPOS (ms):\n");
    printf("  Conversion Listas -> CSR  : %.3f ms\n", t_conv);
    printf("  %d Consultas (Listas)  : %.3f ms\n", consultas, t_cons_listas);
    printf("  %d Consultas (CSR)     : %.3f ms\n", consultas, t_cons_csr);
    printf("  Suma CSR (A + A)          : %.3f ms\n", t_suma);
    printf("  Transpuesta CSR (A^T)     : %.3f ms\n", t_trans);
    printf("  Producto CSR (A * A^T)    : %.3f ms (nnz resultado = %d)\n", t_prod, P ? P->nnz : 0);

    free(qf); free(qc);
    destruirCSR(P); destruirCSR(T); destruirCSR(S); destruirCSR(A);
    liberarListasBench(M);
}

void ejecutarBenchmark(void) {
    int op;
    printf("\n--- MODULO BENCHMARK ---\n");
    printf("[1] Suite oficial (10k, 50k, 100k | 10000x5000 | semilla 42)\n");
    printf("[2] Ingresar parametros personalizados\n");
    printf("Opcion: ");
    if (scanf("%d", &op) != 1) return;

    if (op == 1) {
        const int casos[] = {10000, 50000, 100000};
        for (int i = 0; i < 3; i++)
            ejecutarBenchmarkPersonalizado(10000, 5000, casos[i], 42, 10000);
    } else if (op == 2) {
        int m, n, nnz, consultas;
        unsigned int semilla;
        printf("Filas (m): "); scanf("%d", &m);
        printf("Columnas (n): "); scanf("%d", &n);
        printf("Elementos no nulos (nnz): "); scanf("%d", &nnz);
        printf("Semilla: "); scanf("%u", &semilla);
        printf("Numero de consultas: "); scanf("%d", &consultas);
        ejecutarBenchmarkPersonalizado(m, n, nnz, semilla, consultas);
    }
}