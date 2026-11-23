#ifndef COMPARADOR_OPCIONAL_H
#define COMPARADOR_OPCIONAL_H

#include "matriz_listas.h"
#include "matriz_csr.h"

MatrizDispersa* sumarListas(const MatrizDispersa* A, const MatrizDispersa* B);
MatrizDispersa* transpuestaListas(const MatrizDispersa* A);
void ejecutarComparacionOpcional(void);

#endif