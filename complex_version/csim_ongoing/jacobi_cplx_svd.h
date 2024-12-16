#pragma once 

#include "mat_utils_cplx.h"
#include "jacobi_cplx_evd.h"


void jacobi_svd(Complex* A, Complex* u, Complex* v, Complex* s);

void normalize_matrix(Complex *A, int rows, int cols, T &max_val);
void restore_sigma(Complex *SIGMA, int rows, int cols, T max_val);
