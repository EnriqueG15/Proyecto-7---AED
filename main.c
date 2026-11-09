 /* main.c  -  Integrante 3 (con apoyo del equipo)
 * Menu interactivo que integra listas cruzadas (I1), CSR (I2) y el recomendador (I3).
 * Compilar:
 *   gcc -Wall -O2 main.c matriz_listas.c matriz_csr.c recomendador.c -o proyecto7 -lm
 */
#ifndef _WIN32
#define _POSIX_C_SOURCE 200112L   /* clock_gettime */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#endif

#include "recomendador.h"   /* incluye matriz_listas.h y matriz_csr.h */

/* ====================== Temporizador de alta precision ====================== */

static double ahoraMs(void) {
#ifdef _WIN32
    LARGE_INTEGER f, c;
    QueryPerformanceFrequency(&f);
    QueryPerformanceCounter(&c);
    return 1000.0 * (double)c.QuadPart / (double)f.QuadPart;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
#endif
}

/* ====================== Utilidades de consola ====================== */

static void leerLinea(char* buf, int tam) {
    if (!fgets(buf, tam, stdin)) { printf("\n"); exit(0); }   /* fin de entrada */
}

static int leerEntero(const char* msg, int lo, int hi) {
    char buf[128];
    for (;;) {
        printf("%s", msg);
        fflush(stdout);
        leerLinea(buf, (int)sizeof(buf));
        char* fin;
        long v = strtol(buf, &fin, 10);
        if (fin != buf && (*fin == '\0' || isspace((unsigned char)*fin)) && v >= lo && v <= hi)
            return (int)v;
        printf("  Valor invalido (rango %d..%d).\n", lo, hi);
    }
}

static float leerFloat(const char* msg, float lo, float hi) {
    char buf[128];
    for (;;) {
        printf("%s", msg);
        fflush(stdout);
        leerLinea(buf, (int)sizeof(buf));
        char* fin;
        double v = strtod(buf, &fin);
        if (fin != buf && (*fin == '\0' || isspace((unsigned char)*fin)) && v >= lo && v <= hi)
            return (float)v;
        printf("  Valor invalido (rango %.1f..%.1f).\n", lo, hi);
    }
}

static void leerTexto(const char* msg, const char* porDefecto, char* out, int tam) {
    printf("%s", msg);
    fflush(stdout);
    leerLinea(out, tam);
    out[strcspn(out, "\r\n")] = '\0';
    if (out[0] == '\0') { strncpy(out, porDefecto, (size_t)tam - 1); out[tam - 1] = '\0'; }
}

/* ====================== Estado del sistema ======================
 * L (listas cruzadas) es la fuente de verdad: admite insertar/eliminar.
 * csr y normas son una "vista de calculo" que se reconstruye si hubo cambios. */
typedef struct {
    MatrizDispersa* L;
    MatrizCSR*      csr;
    float*          normas;
    bool            sucio;
} Sistema;

static bool hayDatos(const Sistema* S) {
    if (!S->L) { printf("  Primero carga datos (opcion 1).\n"); return false; }
    return true;
}

static bool refrescar(Sistema* S) {
    if (!hayDatos(S)) return false;
    if (!S->sucio && S->csr && S->normas) return true;
    if (S->csr)    { destruirCSR(S->csr); S->csr = NULL; }
    if (S->normas) { free(S->normas); S->normas = NULL; }
    double t0 = ahoraMs();
    S->csr = convertirListasACSR(S->L);
    if (!S->csr) { printf("  Error: no se pudo convertir a CSR.\n"); return false; }
    S->normas = calcularNormas(S->csr);
    if (!S->normas) { printf("  Error: sin memoria para las normas.\n"); return false; }
    S->sucio = false;
    printf("  [Listas -> CSR + normas: %.3f ms]\n", ahoraMs() - t0);
    return true;
}

static void liberarSistema(Sistema* S) {
    if (S->csr)    destruirCSR(S->csr);
    if (S->normas) free(S->normas);
    if (S->L)      destruirMatrizDispersa(S->L);
    S->csr = NULL; S->normas = NULL; S->L = NULL; S->sucio = true;
}

/* ====================== Opcion 1: cargar datos ====================== */

static void opCargar(Sistema* S) {
    printf("\n  1) MovieLens (u.data)\n  2) Generador sintetico\n  0) Volver\n");
    int o = leerEntero("  Opcion: ", 0, 2);
    if (o == 0) return;

    int m = leerEntero(o == 1 ? "  Usuarios (943 en u.data): " : "  Usuarios m: ", 1, 1000000);
    int n = leerEntero(o == 1 ? "  Items (1682 en u.data): "   : "  Items n: ", 1, 1000000);

    MatrizDispersa* nueva = crearMatriz(m, n);
    if (!nueva) { printf("  Error al crear la matriz.\n"); return; }

    double t0;
    if (o == 1) {
        char ruta[256];
        leerTexto("  Ruta del archivo [u.data]: ", "u.data", ruta, (int)sizeof(ruta));
        t0 = ahoraMs();
        int leidas = cargarMovieLens(nueva, ruta);
        if (leidas < 0) {
            printf("  No se pudo abrir '%s'.\n", ruta);
            destruirMatrizDispersa(nueva);
            return;
        }
        printf("  Calificaciones leidas: %d\n", leidas);
    } else {
        int objetivo = leerEntero("  nnz objetivo (ej. 100000): ", 1, 100000000);
        unsigned int semilla = (unsigned int)leerEntero("  Semilla (0 = por defecto): ", 0, 2000000000);
        t0 = ahoraMs();
        generarSintetico(nueva, objetivo, semilla);
    }
    double ms = ahoraMs() - t0;

    liberarSistema(S);
    S->L = nueva;
    S->sucio = true;

    double celdas = (double)m * (double)n;
    printf("  Matriz %d x %d | nnz = %d | vacias = %.2f %% | carga: %.3f ms\n",
           m, n, S->L->nnz, (1.0 - S->L->nnz / celdas) * 100.0, ms);
    refrescar(S);
}

/* ====================== Opcion 2: CRUD en vivo ====================== */

static void opCRUD(Sistema* S) {
    if (!hayDatos(S)) return;
    printf("\n  1) Consultar  2) Insertar/Actualizar  3) Eliminar  0) Volver\n");
    int o = leerEntero("  Opcion: ", 0, 3);
    if (o == 0) return;

    int u = leerEntero("  Usuario (1-based): ", 1, S->L->m) - 1;
    int i = leerEntero("  Item (1-based): ", 1, S->L->n) - 1;

    if (o == 1) {
        float r = consultarCalificacion(S->L, u, i);
        if (r == 0.0f) printf("  El usuario %d no califico el item %d.\n", u + 1, i + 1);
        else printf("  Calificacion (%d, %d) = %.1f\n", u + 1, i + 1, r);
    } else if (o == 2) {
        float r = leerFloat("  Calificacion (1..5): ", 1.0f, 5.0f);
        float previo = consultarCalificacion(S->L, u, i);
        insertarCalificacion(S->L, u, i, r);
        S->sucio = true;
        printf("  %s nnz = %d\n", previo == 0.0f ? "Insertada." : "Actualizada.", S->L->nnz);
    } else {
        if (consultarCalificacion(S->L, u, i) == 0.0f) { printf("  No existe esa calificacion.\n"); return; }
        eliminarCalificacion(S->L, u, i);
        S->sucio = true;
        printf("  Eliminada. nnz = %d\n", S->L->nnz);
    }
}

/* ====================== Opcion 3: similitud ====================== */

static void opSimilitud(Sistema* S) {
    if (!refrescar(S)) return;
    int u = leerEntero("  Usuario 1 (1-based): ", 1, S->csr->m) - 1;
    int v = leerEntero("  Usuario 2 (1-based): ", 1, S->csr->m) - 1;
    double t0 = ahoraMs();
    float s = similitudCoseno(S->csr, S->normas, u, v);
    double us = (ahoraMs() - t0) * 1000.0;
    printf("  Calificaciones: usuario %d = %d, usuario %d = %d\n",
           u + 1, S->csr->ptr_fila[u + 1] - S->csr->ptr_fila[u],
           v + 1, S->csr->ptr_fila[v + 1] - S->csr->ptr_fila[v]);
    printf("  sim(%d, %d) = %.4f   [%.2f us]\n", u + 1, v + 1, s, us);
}

/* ====================== Opcion 4: Top-k ====================== */

static void opTopK(Sistema* S) {
    if (!refrescar(S)) return;
    int u = leerEntero("  Usuario objetivo (1-based): ", 1, S->csr->m) - 1;
    int k = leerEntero("  k (cuantos items recomendar, 1..100): ", 1, 100);
    int nv = leerEntero("  Vecinos a considerar (1..500, sugerido 30): ", 1, 500);

    ParHeap* vec = (ParHeap*)malloc((size_t)nv * sizeof(ParHeap));
    Recomendacion* rec = (Recomendacion*)malloc((size_t)k * sizeof(Recomendacion));
    if (!vec || !rec) { printf("  Sin memoria.\n"); free(vec); free(rec); return; }

    double t0 = ahoraMs();
    int cv = usuariosSimilares(S->csr, S->normas, u, nv, vec);
    int cr = recomendarTopK(S->csr, S->normas, u, k, nv, rec);
    double ms = ahoraMs() - t0;

    printf("\n  Usuario %d (%d calificaciones)\n", u + 1,
           S->csr->ptr_fila[u + 1] - S->csr->ptr_fila[u]);
    if (cv == 0) {
        printf("  Sin vecinos similares (el usuario no tiene calificaciones o no comparte items).\n");
    } else {
        printf("  Vecinos mas similares:\n");
        for (int j = 0; j < cv && j < 5; j++)
            printf("    usuario %d  sim = %.4f\n", vec[j].id + 1, vec[j].puntaje);
        printf("  Top-%d recomendaciones (%d encontradas):\n", k, cr);
        for (int j = 0; j < cr; j++)
            printf("    %d. item %d  puntaje estimado = %.2f\n", j + 1, rec[j].id_item + 1, rec[j].puntaje);
    }
    printf("  [Tiempo: %.3f ms]\n", ms);
    free(vec); free(rec);
}

/* ====================== Opcion 5: suma / transpuesta / producto ====================== */

static void opAlgebra(Sistema* S) {
    if (!refrescar(S)) return;
    printf("\n  1) Transpuesta A^T\n  2) Suma A + A\n"
           "  3) Producto A x A^T (co-ocurrencias usuario-usuario)\n  0) Volver\n");
    int o = leerEntero("  Opcion: ", 0, 3);
    if (o == 0) return;
    MatrizCSR* A = S->csr;
    int paso = A->nnz > 2000 ? A->nnz / 2000 : 1;   /* se verifican ~2000 entradas */

    if (o == 1) {
        double t0 = ahoraMs();
        MatrizCSR* T = transpuestaCSR(A);
        double ms = ahoraMs() - t0;
        if (!T) { printf("  Error.\n"); return; }
        int fallos = 0, vistos = 0, k = 0;
        for (int u = 0; u < A->m; u++)
            for (int p = A->ptr_fila[u]; p < A->ptr_fila[u + 1]; p++, k++)
                if (k % paso == 0) {
                    vistos++;
                    if (consultarCSR(T, A->columnas[p], u) != A->valores[p]) fallos++;
                }
        printf("  A^T: %d x %d, nnz = %d, %.3f ms\n", T->m, T->n, T->nnz, ms);
        printf("  Verificacion A[u][i] == A^T[i][u] en %d entradas: %s\n", vistos,
               (fallos == 0 && T->nnz == A->nnz) ? "OK" : "FALLO");
        destruirCSR(T);
    } else if (o == 2) {
        double t0 = ahoraMs();
        MatrizCSR* C = sumarCSR(A, A);
        double ms = ahoraMs() - t0;
        if (!C) { printf("  Error.\n"); return; }
        int fallos = 0, vistos = 0, k = 0;
        for (int u = 0; u < A->m; u++)
            for (int p = A->ptr_fila[u]; p < A->ptr_fila[u + 1]; p++, k++)
                if (k % paso == 0) {
                    vistos++;
                    if (fabsf(consultarCSR(C, u, A->columnas[p]) - 2.0f * A->valores[p]) > 1e-4f) fallos++;
                }
        printf("  A + A: %d x %d, nnz = %d, %.3f ms\n", C->m, C->n, C->nnz, ms);
        printf("  Verificacion (A+A)[u][i] == 2*A[u][i] en %d entradas: %s\n", vistos,
               (fallos == 0 && C->nnz == A->nnz) ? "OK" : "FALLO");
        destruirCSR(C);
    } else {
        if (A->m > 2000) {
            printf("  Aviso: A x A^T sera de %d x %d y puede ocupar mucha RAM.\n", A->m, A->m);
            if (leerEntero("  Continuar? (1 = si, 0 = no): ", 0, 1) == 0) return;
        }
        MatrizCSR* T = transpuestaCSR(A);
        if (!T) { printf("  Error al transponer.\n"); return; }
        double t0 = ahoraMs();
        MatrizCSR* C = multiplicarCSR(A, T);
        double ms = ahoraMs() - t0;
        if (!C) { printf("  Error (o sin memoria).\n"); destruirCSR(T); return; }
        printf("  A x A^T: %d x %d, nnz = %d, %.3f ms\n", C->m, C->n, C->nnz, ms);
        int u = leerEntero("  Usuario para comprobar C[u][u] = ||Ru||^2 (1-based): ", 1, A->m) - 1;
        float diag = consultarCSR(C, u, u);
        float esperado = S->normas[u] * S->normas[u];
        printf("  C[u][u] = %.2f  |  ||Ru||^2 = %.2f  -> %s\n", diag, esperado,
               fabsf(diag - esperado) <= 1e-3f * (1.0f + esperado) ? "OK" : "FALLO");
        destruirCSR(C);
        destruirCSR(T);
    }
}

/* ====================== Opcion 6: memoria y tiempos ====================== */

static void opReporte(Sistema* S) {
    if (!refrescar(S)) return;
    MatrizCSR* A = S->csr;
    size_t m = (size_t)A->m, n = (size_t)A->n, nnz = (size_t)A->nnz;

    /* Reporte exacto del Integrante 1 (listas vs densa) */
    imprimirReporteMemoria(S->L);

    /* Comparativa de las tres representaciones */
    size_t bDensa  = m * n * sizeof(float) + m * sizeof(float*);
    size_t bListas = nnz * sizeof(Nodo) + (m + n) * sizeof(Nodo*) + sizeof(MatrizDispersa);
    size_t bCSR    = memoriaUsadaCSR(A);   /* funcion del Integrante 2 */
    printf("  COMPARATIVA DE MEMORIA (bytes)\n");
    printf("    Densa  : %zu (%.2f MB)\n", bDensa,  bDensa  / (1024.0 * 1024.0));
    printf("    Listas : %zu (%.2f MB)\n", bListas, bListas / (1024.0 * 1024.0));
    printf("    CSR    : %zu (%.2f MB)\n", bCSR,    bCSR    / (1024.0 * 1024.0));

    /* Tiempos */
    printf("\n  TIEMPOS DE EJECUCION\n");
    unsigned int s = 12345u;
    const int NQ = 100000;
    double acc = 0.0;

    double t0 = ahoraMs();
    for (int q = 0; q < NQ; q++) {
        s ^= s << 13; s ^= s >> 17; s ^= s << 5;
        int u = (int)(s % (unsigned)A->m);
        s ^= s << 13; s ^= s >> 17; s ^= s << 5;
        int i = (int)(s % (unsigned)A->n);
        acc += consultarCalificacion(S->L, u, i);
    }
    double msListas = ahoraMs() - t0;

    s = 12345u;
    t0 = ahoraMs();
    for (int q = 0; q < NQ; q++) {
        s ^= s << 13; s ^= s >> 17; s ^= s << 5;
        int u = (int)(s % (unsigned)A->m);
        s ^= s << 13; s ^= s >> 17; s ^= s << 5;
        int i = (int)(s % (unsigned)A->n);
        acc += consultarCSR(A, u, i);
    }
    double msCSR = ahoraMs() - t0;
    printf("    %d consultas aleatorias: Listas = %.3f ms | CSR = %.3f ms (checksum %.0f)\n",
           NQ, msListas, msCSR, acc);

    t0 = ahoraMs();
    MatrizCSR* T = transpuestaCSR(A);
    printf("    Transpuesta: %.3f ms\n", ahoraMs() - t0);
    if (T) destruirCSR(T);

    const int NS = 2000;
    double sumaSim = 0.0;
    t0 = ahoraMs();
    for (int q = 0; q < NS; q++) {
        s ^= s << 13; s ^= s >> 17; s ^= s << 5;
        int u = (int)(s % (unsigned)A->m);
        s ^= s << 13; s ^= s >> 17; s ^= s << 5;
        int v = (int)(s % (unsigned)A->m);
        sumaSim += similitudCoseno(A, S->normas, u, v);
    }
    printf("    %d similitudes coseno: %.3f ms (sim media %.4f)\n", NS, ahoraMs() - t0, sumaSim / NS);

    Recomendacion rec[10];
    t0 = ahoraMs();
    for (int q = 0; q < 5; q++) recomendarTopK(A, S->normas, (q * 7) % A->m, 10, 30, rec);
    printf("    Top-10 (30 vecinos), promedio de 5 usuarios: %.3f ms\n\n", (ahoraMs() - t0) / 5.0);
}

/* ====================== Menu principal ====================== */

int main(void) {
    Sistema S = { NULL, NULL, NULL, true };

    for (;;) {
        printf("\n=========== PROYECTO 7: RECOMENDADOR CON MATRICES DISPERSAS ===========\n");
        if (S.L) printf("  [Matriz %d x %d, nnz = %d]\n", S.L->m, S.L->n, S.L->nnz);
        else     printf("  [Sin datos cargados]\n");
        printf("  1) Cargar calificaciones (MovieLens / sintetico)\n"
               "  2) Consultar / Insertar / Eliminar una calificacion\n"
               "  3) Similitud coseno entre dos usuarios\n"
               "  4) Top-k recomendaciones para un usuario\n"
               "  5) Suma, transpuesta y producto\n"
               "  6) Informe de memoria (Listas vs CSR vs Densa) y tiempos\n"
               "  0) Salir\n");
        int op = leerEntero("  Opcion: ", 0, 6);
        switch (op) {
            case 1: opCargar(&S);     break;
            case 2: opCRUD(&S);       break;
            case 3: opSimilitud(&S);  break;
            case 4: opTopK(&S);       break;
            case 5: opAlgebra(&S);    break;
            case 6: opReporte(&S);    break;
            case 0: liberarSistema(&S); printf("  Hasta luego.\n"); return 0;
        }
    }
}