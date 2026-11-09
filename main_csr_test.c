/*
 * Pruebas independientes del modulo CSR (Integrante 2).
 * No usa funciones del Integrante 1: construye las listas cruzadas de prueba
 * con funciones locales (static) sobre las estructuras de matriz_listas.h.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "matriz_listas.h"
#include "matriz_csr.h"

#define TOL 1e-4f

static int g_total = 0, g_ok = 0;

static float absf(float x) { return x < 0.0f ? -x : x; }

static void verificar(const char* descripcion, bool condicion) {
    g_total++;
    if (condicion) g_ok++;
    printf("  [%s] %s\n", condicion ? " OK  " : "FALLO", descripcion);
}

/* =============== Listas cruzadas de prueba (solo para tests) =============== */

static MatrizDispersa* crearListasPrueba(int m, int n) {
    MatrizDispersa* M = (MatrizDispersa*)malloc(sizeof(MatrizDispersa));
    if (M == NULL) return NULL;
    M->m = m;
    M->n = n;
    M->nnz = 0;
    M->cabeceras_fila = (Nodo**)calloc((size_t)m, sizeof(Nodo*));
    M->cabeceras_col  = (Nodo**)calloc((size_t)n, sizeof(Nodo*));
    M->avail = NULL;
    if (M->cabeceras_fila == NULL || M->cabeceras_col == NULL) {
        free(M->cabeceras_fila); free(M->cabeceras_col); free(M);
        return NULL;
    }
    return M;
}

/* Insercion ordenada en fila (sig_col) y columna (sig_fila). Sobrescribe si existe. */
static bool insertarListasPrueba(MatrizDispersa* M, int i, int j, float v) {
    if (M == NULL || i < 0 || i >= M->m || j < 0 || j >= M->n || v == 0.0f) return false;

    Nodo* ant = NULL;
    Nodo* p = M->cabeceras_fila[i];
    while (p != NULL && p->col < j) { ant = p; p = p->sig_col; }
    if (p != NULL && p->col == j) { p->calificacion = v; return true; }

    Nodo* q = (Nodo*)malloc(sizeof(Nodo));
    if (q == NULL) return false;
    q->fila = i; q->col = j; q->calificacion = v;
    q->sig_col = p;
    if (ant != NULL) ant->sig_col = q; else M->cabeceras_fila[i] = q;

    Nodo* antc = NULL;
    Nodo* r = M->cabeceras_col[j];
    while (r != NULL && r->fila < i) { antc = r; r = r->sig_fila; }
    q->sig_fila = r;
    if (antc != NULL) antc->sig_fila = q; else M->cabeceras_col[j] = q;

    M->nnz++;
    return true;
}

static void liberarListasPrueba(MatrizDispersa* M) {
    if (M == NULL) return;
    for (int i = 0; i < M->m; i++) {
        Nodo* p = M->cabeceras_fila[i];
        while (p != NULL) { Nodo* s = p->sig_col; free(p); p = s; }
    }
    Nodo* p = M->avail;
    while (p != NULL) { Nodo* s = p->sig_col; free(p); p = s; }
    free(M->cabeceras_fila);
    free(M->cabeceras_col);
    free(M);
}

static MatrizDispersa* construirDesdeTripletas(int m, int n, const int* fi, const int* co,
                                               const float* va, int k) {
    MatrizDispersa* M = crearListasPrueba(m, n);
    if (M == NULL) return NULL;
    for (int t = 0; t < k; t++) insertarListasPrueba(M, fi[t], co[t], va[t]);
    return M;
}

/* ===================== Referencias densas (oraculo) ===================== */

static float* listasADensa(const MatrizDispersa* M) {
    float* D = (float*)calloc((size_t)M->m * M->n, sizeof(float));
    if (D == NULL) return NULL;
    for (int i = 0; i < M->m; i++)
        for (const Nodo* p = M->cabeceras_fila[i]; p != NULL; p = p->sig_col)
            D[(size_t)i * M->n + p->col] = p->calificacion;
    return D;
}

static float* csrADensa(const MatrizCSR* A) {
    float* D = (float*)calloc((size_t)A->m * A->n, sizeof(float));
    if (D == NULL) return NULL;
    for (int i = 0; i < A->m; i++)
        for (int k = A->ptr_fila[i]; k < A->ptr_fila[i + 1]; k++)
            D[(size_t)i * A->n + A->columnas[k]] = A->valores[k];
    return D;
}

static float* sumaDensa(const float* X, const float* Y, int m, int n) {
    float* R = (float*)calloc((size_t)m * n, sizeof(float));
    if (R == NULL) return NULL;
    for (int t = 0; t < m * n; t++) R[t] = X[t] + Y[t];
    return R;
}

static float* transpuestaDensa(const float* X, int m, int n) {
    float* R = (float*)calloc((size_t)m * n, sizeof(float));
    if (R == NULL) return NULL;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            R[(size_t)j * m + i] = X[(size_t)i * n + j];
    return R;
}

static float* productoDenso(const float* X, const float* Y, int m, int k, int n) {
    float* R = (float*)calloc((size_t)m * n, sizeof(float));
    if (R == NULL) return NULL;
    for (int i = 0; i < m; i++)
        for (int t = 0; t < k; t++) {
            float a = X[(size_t)i * k + t];
            if (a == 0.0f) continue;
            for (int j = 0; j < n; j++) R[(size_t)i * n + j] += a * Y[(size_t)t * n + j];
        }
    return R;
}

/* Compara un CSR con una referencia densa de dimensiones m x n */
static bool coincideConDensa(const MatrizCSR* A, const float* D, int m, int n) {
    if (A == NULL || D == NULL || A->m != m || A->n != n) return false;
    float* X = csrADensa(A);
    if (X == NULL) return false;
    bool ok = true;
    for (int t = 0; t < m * n && ok; t++)
        if (absf(X[t] - D[t]) > TOL) ok = false;
    free(X);
    return ok;
}

/* Invariantes estructurales del formato CSR */
static bool invariantesCSR(const MatrizCSR* A) {
    if (A == NULL) return false;
    if (A->ptr_fila[0] != 0 || A->ptr_fila[A->m] != A->nnz) return false;
    for (int i = 0; i < A->m; i++) {
        if (A->ptr_fila[i] > A->ptr_fila[i + 1]) return false;
        for (int k = A->ptr_fila[i]; k < A->ptr_fila[i + 1]; k++) {
            if (A->columnas[k] < 0 || A->columnas[k] >= A->n) return false;
            if (k > A->ptr_fila[i] && A->columnas[k] <= A->columnas[k - 1]) return false;
            if (A->valores[k] == 0.0f) return false;     /* no se guardan ceros */
        }
    }
    return true;
}

/* ============================= Presentacion ============================= */

static void mostrarCSR(const char* nombre, const MatrizCSR* A) {
    if (A == NULL) { printf("\n--- %s: (no disponible) ---\n", nombre); return; }
    printf("\n--- %s  (%d x %d, nnz = %d, memoria = %lu bytes) ---\n",
           nombre, A->m, A->n, A->nnz, (unsigned long)memoriaUsadaCSR(A));
    printf("Estructura interna CSR:\n");
    imprimirEstructuraCSR(A);
    printf("Tripletas (fila, columna, valor):\n");
    imprimirTripletasCSR(A);
    printf("Vista matricial:\n");
    imprimirVistaMatricialCSR(A);
}

/* ========================= Pruebas automaticas ========================= */

static void ejecutarPruebasAutomaticas(void) {
    g_total = 0;
    g_ok = 0;

    /* A (4x5): matriz usuarios x items. La fila 2 esta vacia a proposito. */
    const int   aF[] = {0, 0, 1, 1, 3, 3, 3};
    const int   aC[] = {0, 3, 1, 4, 0, 2, 4};
    const float aV[] = {5, 3, 4, 1, 2, 5, 4};

    /* B (4x5): (0,0) y (3,4) cancelan con A para probar la omision de ceros */
    const int   bF[] = {0, 0, 1, 2, 3};
    const int   bC[] = {0, 2, 1, 3, 4};
    const float bV[] = {-5, 1, 1, 2, -4};

    /* P (5x3): para el producto A(4x5) * P(5x3) */
    const int   pF[] = {0, 0, 1, 2, 3, 4};
    const int   pC[] = {0, 2, 1, 0, 2, 1};
    const float pV[] = {1, 2, 3, 1, 4, 2};

    MatrizDispersa* LA = construirDesdeTripletas(4, 5, aF, aC, aV, 7);
    MatrizDispersa* LB = construirDesdeTripletas(4, 5, bF, bC, bV, 5);
    MatrizDispersa* LP = construirDesdeTripletas(5, 3, pF, pC, pV, 6);
    MatrizDispersa* LV = crearListasPrueba(3, 3);   /* matriz vacia */
    if (!LA || !LB || !LP || !LV) {
        printf("[ERROR] Memoria insuficiente al construir las listas de prueba.\n");
        liberarListasPrueba(LA); liberarListasPrueba(LB);
        liberarListasPrueba(LP); liberarListasPrueba(LV);
        return;
    }

    float* DA = listasADensa(LA);
    float* DB = listasADensa(LB);
    float* DP = listasADensa(LP);

    printf("\n================ PRUEBA 1: Conversion Listas -> CSR ================\n");
    MatrizCSR* A = convertirListasACSR(LA);
    MatrizCSR* B = convertirListasACSR(LB);
    MatrizCSR* P = convertirListasACSR(LP);
    MatrizCSR* V = convertirListasACSR(LV);
    mostrarCSR("A", A);
    mostrarCSR("B", B);
    mostrarCSR("P", P);
    verificar("A convertida coincide con las listas cruzadas", coincideConDensa(A, DA, 4, 5));
    verificar("B convertida coincide con las listas cruzadas", coincideConDensa(B, DB, 4, 5));
    verificar("P convertida coincide con las listas cruzadas", coincideConDensa(P, DP, 5, 3));
    verificar("nnz(A) = 7 y ptr_fila[4] = 7", A && A->nnz == 7 && A->ptr_fila[4] == 7);
    verificar("Fila vacia 2: ptr_fila[2] == ptr_fila[3]", A && A->ptr_fila[2] == A->ptr_fila[3]);
    verificar("Invariantes CSR de A, B y P", invariantesCSR(A) && invariantesCSR(B) && invariantesCSR(P));
    verificar("Matriz vacia: nnz = 0 y ptr_fila todo en 0",
              V && V->nnz == 0 && V->ptr_fila[0] == 0 && V->ptr_fila[3] == 0);
    verificar("Conversion de NULL retorna NULL", convertirListasACSR(NULL) == NULL);

    printf("\n================ PRUEBA 2: Consulta por busqueda binaria ================\n");
    bool todasOk = true;
    for (int i = 0; i < 4 && A; i++)
        for (int j = 0; j < 5; j++)
            if (absf(consultarCSR(A, i, j) - DA[i * 5 + j]) > TOL) todasOk = false;
    printf("  consultarCSR(A, 0, 3) = %.2f (esperado 3.00)\n", consultarCSR(A, 0, 3));
    printf("  consultarCSR(A, 3, 4) = %.2f (esperado 4.00)\n", consultarCSR(A, 3, 4));
    printf("  consultarCSR(A, 2, 2) = %.2f (esperado 0.00, fila vacia)\n", consultarCSR(A, 2, 2));
    printf("  consultarCSR(A, 1, 2) = %.2f (esperado 0.00, hueco)\n", consultarCSR(A, 1, 2));
    verificar("Las 20 posiciones de A coinciden con la referencia densa", todasOk);
    verificar("Indices fuera de rango retornan 0",
              consultarCSR(A, -1, 0) == 0.0f && consultarCSR(A, 4, 0) == 0.0f &&
              consultarCSR(A, 0, 5) == 0.0f && consultarCSR(NULL, 0, 0) == 0.0f);

    printf("\n================ PRUEBA 3: Suma A + B (merge de 2 punteros) ================\n");
    MatrizCSR* S = sumarCSR(A, B);
    float* DS = sumaDensa(DA, DB, 4, 5);
    mostrarCSR("A + B", S);
    verificar("A + B coincide con la suma densa", coincideConDensa(S, DS, 4, 5));
    verificar("Cancelaciones (0,0) y (3,4) no se almacenan (nnz = 7)", S && S->nnz == 7);
    verificar("Invariantes CSR de A + B", invariantesCSR(S));
    MatrizCSR* AV = sumarCSR(A, P);
    verificar("Suma con dimensiones distintas retorna NULL", AV == NULL);
    destruirCSR(AV);

    printf("\n================ PRUEBA 4: Transpuesta A^T (Counting Sort) ================\n");
    MatrizCSR* T  = transpuestaCSR(A);
    MatrizCSR* TT = transpuestaCSR(T);
    float* DT = transpuestaDensa(DA, 4, 5);
    mostrarCSR("A^T", T);
    verificar("A^T coincide con la transpuesta densa (5 x 4)", coincideConDensa(T, DT, 5, 4));
    verificar("Invariantes CSR de A^T (columnas ordenadas)", invariantesCSR(T));
    verificar("(A^T)^T es identica a A", sonIgualesCSR(TT, A, TOL));
    MatrizCSR* VT = transpuestaCSR(V);
    verificar("Transpuesta de matriz vacia: nnz = 0", VT && VT->nnz == 0);
    destruirCSR(VT);

    printf("\n================ PRUEBA 5: Producto (Gustavson) ================\n");
    MatrizCSR* AP = multiplicarCSR(A, P);
    float* DAP = productoDenso(DA, DP, 4, 5, 3);
    mostrarCSR("A * P", AP);
    verificar("A * P (4x3) coincide con el producto denso", coincideConDensa(AP, DAP, 4, 3));
    verificar("Invariantes CSR de A * P", invariantesCSR(AP));

    MatrizCSR* AAT = multiplicarCSR(A, T);
    float* DAAT = productoDenso(DA, DT, 4, 5, 4);
    mostrarCSR("A * A^T (similitud usuario-usuario)", AAT);
    verificar("A * A^T (4x4) coincide con el producto denso", coincideConDensa(AAT, DAAT, 4, 4));
    MatrizCSR* AATt = transpuestaCSR(AAT);
    verificar("A * A^T == (A * A^T)^T", sonIgualesCSR(AAT, AATt, TOL));
    destruirCSR(AATt);

    MatrizCSR* AB = multiplicarCSR(A, B);
    verificar("Producto incompatible (4x5 * 4x5) retorna NULL", AB == NULL);
    destruirCSR(AB);

    printf("\n================ PRUEBA 6: Memoria ================\n");
    if (A != NULL) {
        size_t esperado = sizeof(MatrizCSR) + 7 * sizeof(float) + 7 * sizeof(int) + 5 * sizeof(int);
        printf("  memoriaUsadaCSR(A) = %lu bytes\n", (unsigned long)memoriaUsadaCSR(A));
        printf("  Densa equivalente   = %lu bytes\n", (unsigned long)(4 * 5 * sizeof(float)));
        verificar("memoriaUsadaCSR(A) = sizeof + 7*float + 7*int + 5*int",
                  memoriaUsadaCSR(A) == esperado);
    }
    verificar("memoriaUsadaCSR(NULL) = 0", memoriaUsadaCSR(NULL) == 0);

    printf("\n=====================================================================\n");
    printf("  RESULTADO: %d / %d pruebas superadas %s\n", g_ok, g_total,
           g_ok == g_total ? "(TODO CORRECTO)" : "(REVISAR FALLOS)");
    printf("=====================================================================\n");

    free(DA); free(DB); free(DP); free(DS); free(DT); free(DAP); free(DAAT);
    destruirCSR(A); destruirCSR(B); destruirCSR(P); destruirCSR(V);
    destruirCSR(S); destruirCSR(T); destruirCSR(TT); destruirCSR(AP); destruirCSR(AAT);
    liberarListasPrueba(LA); liberarListasPrueba(LB);
    liberarListasPrueba(LP); liberarListasPrueba(LV);
}

/* ========================= Modo manual (teclado) ========================= */

static void limpiarBuffer(void) {
    int c;
    while ((c = getchar()) != '\n' && c != EOF) { }
}

static int leerEntero(const char* msg, int min, int max) {
    int x, r;
    for (;;) {
        printf("%s", msg);
        fflush(stdout);
        r = scanf("%d", &x);
        if (r == EOF) { printf("\nFin de entrada.\n"); exit(EXIT_SUCCESS); }
        limpiarBuffer();
        if (r == 1 && x >= min && x <= max) return x;
        printf("  [!] Ingrese un entero en [%d, %d].\n", min, max);
    }
}

static MatrizDispersa* leerMatrizTeclado(const char* nombre) {
    printf("\n--- Matriz %s ---\n", nombre);
    int m = leerEntero("  Filas m (1..1000): ", 1, 1000);
    int n = leerEntero("  Columnas n (1..1000): ", 1, 1000);
    int maxk = m * n < 10000 ? m * n : 10000;
    char msg[96];
    snprintf(msg, sizeof msg, "  Cantidad de datos no nulos (0..%d): ", maxk);
    int k = leerEntero(msg, 0, maxk);

    MatrizDispersa* M = crearListasPrueba(m, n);
    if (M == NULL) { printf("  [ERROR] Memoria insuficiente.\n"); return NULL; }
    if (k > 0) printf("  Ingrese: fila columna valor  (indices desde 0)\n");
    for (int t = 0; t < k; t++) {
        int i, j, r;
        float v;
        printf("  [%d/%d] > ", t + 1, k);
        fflush(stdout);
        r = scanf("%d %d %f", &i, &j, &v);
        if (r == EOF) { liberarListasPrueba(M); printf("\nFin de entrada.\n"); exit(EXIT_SUCCESS); }
        limpiarBuffer();
        if (r != 3 || i < 0 || i >= m || j < 0 || j >= n) {
            printf("  [!] Dato invalido, repita.\n");
            t--;
            continue;
        }
        if (v == 0.0f) printf("  (valor 0 ignorado: no se almacena)\n");
        else insertarListasPrueba(M, i, j, v);
    }
    return M;
}

static void ejecutarModoManual(void) {
    MatrizDispersa* LA = leerMatrizTeclado("A");
    if (LA == NULL) return;
    MatrizCSR* A = convertirListasACSR(LA);
    mostrarCSR("A (convertida a CSR)", A);

    /* Consultas */
    printf("\n  Consultas binarias en A (ingrese -1 para terminar):\n");
    for (;;) {
        int i = leerEntero("  Fila i (-1 = fin): ", -1, A->m - 1);
        if (i < 0) break;
        int j = leerEntero("  Columna j: ", 0, A->n - 1);
        printf("  A[%d][%d] = %.2f\n", i, j, consultarCSR(A, i, j));
    }

    MatrizCSR* T = transpuestaCSR(A);
    mostrarCSR("A^T", T);

    int conB = leerEntero("\n  Ingresar matriz B para suma/producto? (1 = si, 0 = no): ", 0, 1);
    if (conB) {
        MatrizDispersa* LB = leerMatrizTeclado("B");
        if (LB != NULL) {
            MatrizCSR* B = convertirListasACSR(LB);
            mostrarCSR("B (convertida a CSR)", B);

            MatrizCSR* S = sumarCSR(A, B);
            if (S) mostrarCSR("A + B", S);
            else printf("\n  A + B: dimensiones incompatibles (%dx%d vs %dx%d).\n",
                        A->m, A->n, B->m, B->n);

            MatrizCSR* P = multiplicarCSR(A, B);
            if (P) mostrarCSR("A * B", P);
            else printf("\n  A * B: incompatibles (A tiene %d columnas, B tiene %d filas).\n",
                        A->n, B->m);

            destruirCSR(S); destruirCSR(P); destruirCSR(B);
            liberarListasPrueba(LB);
        }
    }
    destruirCSR(T); destruirCSR(A);
    liberarListasPrueba(LA);
}

/* ================================ main ================================ */

int main(void) {
    int op;
    printf("==========================================================\n");
    printf("  Proyecto 7 - Pruebas del modulo CSR (Integrante 2)\n");
    printf("==========================================================\n");
    do {
        printf("\n  [1] Ejecutar bateria de pruebas automaticas\n");
        printf("  [2] Ingresar matrices manualmente y probar operaciones\n");
        printf("  [0] Salir\n");
        op = leerEntero("  Opcion: ", 0, 2);
        if (op == 1) ejecutarPruebasAutomaticas();
        else if (op == 2) ejecutarModoManual();
    } while (op != 0);
    printf("Fin de las pruebas.\n");
    return 0;
}