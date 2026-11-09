#ifndef RECOMENDADOR_H
#define RECOMENDADOR_H

/* Recomendador.h  -  Integrante 3 
 * TAD ColaPrioridad (heap binario en arreglo), lector MovieLens / generador
 * sintetico, similitud coseno y recomendacion Top-k.
 * Convencion del equipo: indices base cero (usuarios 0..m-1, items 0..n-1).
 */

#include <stddef.h>
#include <stdbool.h>
#include "matriz_listas.h"
#include "matriz_csr.h"

/* Prototipos del Integrante 1: su matriz_listas.h solo trae los structs.
 * Cuando el los agregue a su header, este bloque puede borrarse. */
MatrizDispersa* crearMatriz(int m, int n);
void  destruirMatrizDispersa(MatrizDispersa* M);
void  insertarCalificacion(MatrizDispersa* M, int i, int j, float val);
float consultarCalificacion(MatrizDispersa* M, int i, int j);
void  eliminarCalificacion(MatrizDispersa* M, int i, int j);
void  imprimirReporteMemoria(MatrizDispersa* M);

/* ====================== TAD Cola de Prioridad ====================== */

typedef enum { HEAP_MAX = 0, HEAP_MIN = 1 } ModoHeap;

typedef struct {
    int   id;        /* id_item (o id de usuario) */
    float puntaje;   /* puntaje estimado / similitud */
} ParHeap;

typedef struct {
    ParHeap* arreglo;
    int      tamano;
    int      capacidad;
    ModoHeap modo;   /* HEAP_MAX: raiz = mayor puntaje. HEAP_MIN: raiz = menor puntaje */
} TColaPrioridad;

/* Pre: capacidad > 0. Post: cola vacia (tamano = 0). Devuelve NULL si falla. */
TColaPrioridad* crearColaPrioridad(int capacidad, ModoHeap modo);
void            destruirColaPrioridad(TColaPrioridad* C);

bool    esVacia(const TColaPrioridad* C);
/* Pre: tamano < capacidad. Devuelve false (sin modificar) si esta llena. O(log n) */
bool    insertarHeap(TColaPrioridad* C, int id, float puntaje);
/* Pre: !esVacia(C). Extrae la raiz y restaura el orden con hundir. O(log n) */
ParHeap extraerOptimo(TColaPrioridad* C);
/* Pre: !esVacia(C). O(1) */
ParHeap verOptimo(const TColaPrioridad* C);
/* Pre: !esVacia(C). Sustituye la raiz y hunde. O(log n) */
void    reemplazarOptimo(TColaPrioridad* C, int id, float puntaje);

/* ====================== Dataset ====================== */

/* Lee u.data de MovieLens 100k (user item rating timestamp) o CSV equivalente.
 * Resta 1 a los ids. Devuelve cuantas calificaciones leyo, o -1 si no abre el archivo. */
int cargarMovieLens(MatrizDispersa* M, const char* ruta);

/* Rellena M con nnzObjetivo calificaciones aleatorias (1..5). Limita a 5% de
 * celdas ocupadas (>95% vacias). Devuelve el nnz resultante. */
int generarSintetico(MatrizDispersa* M, int nnzObjetivo, unsigned int semilla);

/* ====================== Similitud y recomendacion ====================== */

typedef struct {
    int   id_item;   /* base cero */
    float puntaje;
} Recomendacion;

/* Arreglo de m normas ||Ru||, calculado en O(nnz). El llamador hace free(). */
float* calcularNormas(const MatrizCSR* A);

/* sim(u,v) = (Ru . Rv)/(||Ru|| ||Rv||) por interseccion de dos punteros. O(nnz_u + nnz_v) */
float similitudCoseno(const MatrizCSR* A, const float* normas, int u, int v);

/* Hasta nVecinos usuarios mas parecidos a u (sim > 0), ordenados de mayor a menor.
 * salida debe tener espacio para nVecinos. Devuelve cuantos encontro. */
int usuariosSimilares(const MatrizCSR* A, const float* normas, int u, int nVecinos, ParHeap* salida);

/* Top-k items no calificados por u, con puntaje = promedio ponderado por similitud
 * de las calificaciones de sus vecinos. salida debe tener espacio para k.
 * Devuelve cuantos items recomendo (<= k). */
int recomendarTopK(const MatrizCSR* A, const float* normas, int u, int k,
                   int nVecinos, Recomendacion* salida);

#endif