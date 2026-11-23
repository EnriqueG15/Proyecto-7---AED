#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "comparador_opcional.h"

MatrizDispersa* sumarListas(const MatrizDispersa* A, const MatrizDispersa* B) {
    if (!A || !B || A->m != B->m || A->n != B->n) return NULL;
    MatrizDispersa* C = crearMatriz(A->m, A->n);
    if (!C) return NULL;

    for (int i = 0; i < A->m; i++) {
        const Nodo* pa = A->cabeceras_fila[i];
        const Nodo* pb = B->cabeceras_fila[i];
        while (pa != NULL && pb != NULL) {
            if (pa->col < pb->col) {
                insertarCalificacion(C, i, pa->col, pa->calificacion);
                pa = pa->sig_col;
            } else if (pb->col < pa->col) {
                insertarCalificacion(C, i, pb->col, pb->calificacion);
                pb = pb->sig_col;
            } else {
                float s = pa->calificacion + pb->calificacion;
                if (s != 0.0f) insertarCalificacion(C, i, pa->col, s);
                pa = pa->sig_col;
                pb = pb->sig_col;
            }
        }
        while (pa != NULL) { insertarCalificacion(C, i, pa->col, pa->calificacion); pa = pa->sig_col; }
        while (pb != NULL) { insertarCalificacion(C, i, pb->col, pb->calificacion); pb = pb->sig_col; }
    }
    return C;
}

MatrizDispersa* transpuestaListas(const MatrizDispersa* A) {
    if (!A) return NULL;
    MatrizDispersa* T = crearMatriz(A->n, A->m);
    if (!T) return NULL;

    for (int i = 0; i < A->m; i++) {
        for (const Nodo* p = A->cabeceras_fila[i]; p != NULL; p = p->sig_col) {
            insertarCalificacion(T, p->col, i, p->calificacion);
        }
    }
    return T;
}

static unsigned int r32(void) {
    return ((unsigned int)rand() << 15) ^ (unsigned int)rand();
}

static double msDiff(clock_t a, clock_t b) {
    return (double)(b - a) * 1000.0 / CLOCKS_PER_SEC;
}

static void correrCaso(int m, int n, int nnz, unsigned int semilla, int q_count) {
    srand(semilla);
    MatrizDispersa* L = crearMatriz(m, n);
    if (!L) return;

    while (L->nnz < nnz) {
        int i = (int)(r32() % (unsigned int)m);
        int j = (int)(r32() % (unsigned int)n);
        float v = (float)((rand() % 5) + 1);
        insertarCalificacion(L, i, j, v);
    }

    int* qf = (int*)malloc((size_t)q_count * sizeof(int));
    int* qc = (int*)malloc((size_t)q_count * sizeof(int));
    if (!qf || !qc) { free(qf); free(qc); destruirMatrizDispersa(L); return; }

    for (int k = 0; k < q_count; k++) {
        qf[k] = (int)(r32() % (unsigned int)m);
        qc[k] = (int)(r32() % (unsigned int)n);
    }

    clock_t t0 = clock();
    MatrizCSR* C = convertirListasACSR(L);
    clock_t t1 = clock();
    double t_conv = msDiff(t0, t1);

    volatile float acc = 0.0f;
    t0 = clock();
    for (int k = 0; k < q_count; k++) acc += consultarCalificacion(L, qf[k], qc[k]);
    t1 = clock();
    double t_ql = msDiff(t0, t1);

    t0 = clock();
    for (int k = 0; k < q_count; k++) acc += consultarCSR(C, qf[k], qc[k]);
    t1 = clock();
    double t_qc = msDiff(t0, t1);
    (void)acc;

    t0 = clock();
    MatrizDispersa* TL = transpuestaListas(L);
    t1 = clock();
    double t_tl = msDiff(t0, t1);

    t0 = clock();
    MatrizCSR* TC = transpuestaCSR(C);
    t1 = clock();
    double t_tc = msDiff(t0, t1);

    t0 = clock();
    MatrizDispersa* SL = sumarListas(L, L);
    t1 = clock();
    double t_sl = msDiff(t0, t1);

    t0 = clock();
    MatrizCSR* SC = sumarCSR(C, C);
    t1 = clock();
    double t_sc = msDiff(t0, t1);

    size_t b_listas = sizeof(MatrizDispersa) + (size_t)(m + n) * sizeof(Nodo*) + (size_t)L->nnz * sizeof(Nodo);
    size_t b_csr    = memoriaUsadaCSR(C);

    printf("\n[%d x %d | nnz=%d | semilla=%u]\n", m, n, nnz, semilla);
    printf("  Memoria -> Listas: %.2f MB (%zu B) | CSR: %.2f MB (%zu B)\n",
           (double)b_listas / (1024.0 * 1024.0), b_listas,
           (double)b_csr / (1024.0 * 1024.0), b_csr);
    printf("  Conversion Listas->CSR : %.3f ms\n", t_conv);
    printf("  %d consultas           : Listas %.3f ms | CSR %.3f ms\n", q_count, t_ql, t_qc);
    printf("  Transpuesta            : Listas %.3f ms | CSR %.3f ms\n", t_tl, t_tc);
    printf("  Suma (A+A)             : Listas %.3f ms | CSR %.3f ms\n", t_sl, t_sc);

    destruirCSR(SC); destruirCSR(TC); destruirCSR(C);
    destruirMatrizDispersa(SL); destruirMatrizDispersa(TL); destruirMatrizDispersa(L);
    free(qf); free(qc);
}

void ejecutarComparacionOpcional(void) {
    int op;
    printf("\n=== COMPARACION OPCIONAL: LISTAS VS CSR ===\n");
    printf("[1] Suite 10k, 50k, 100k (10000x5000, semilla 42)\n");
    printf("[2] Ingresar datos personalizados\n");
    printf("Opcion: ");
    if (scanf("%d", &op) != 1) return;

    if (op == 1) {
        int casos[3] = {10000, 50000, 100000};
        for (int i = 0; i < 3; i++) correrCaso(10000, 5000, casos[i], 42, 10000);
    } else if (op == 2) {
        int m, n, nnz, q;
        unsigned int sem;
        printf("Filas (m): "); if (scanf("%d", &m) != 1) return;
        printf("Columnas (n): "); if (scanf("%d", &n) != 1) return;
        printf("nnz: "); if (scanf("%d", &nnz) != 1) return;
        printf("Semilla: "); if (scanf("%u", &sem) != 1) return;
        printf("Consultas: "); if (scanf("%d", &q) != 1) return;
        if (m > 0 && n > 0 && nnz >= 0 && (long long)nnz <= (long long)m * n && q >= 0)
            correrCaso(m, n, nnz, sem, q);
        else
            printf("Datos fuera de rango.\n");
    }
}