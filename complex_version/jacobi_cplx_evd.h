#pragma once 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "jacobi_type.h"
#include "jacobi_config.h"
#include "mat_utils_cplx.h"

double complex_norm(Complex z); 
void print_matrix(Complex*A, int N);
void make_identity(Complex* u, int N);
void transpose(Complex* A, Complex* A_t, int N);
void matmul(Complex* A, Complex* B, Complex* C, int N);
void jacobi_rotate(Complex *A, Complex *U, int N, int k, int l);
void find_maxDiagOff(Complex* A, int N, int* p, int *q, double* value);
double cal_sumDiagOff(Complex *A, int N);

void jacobi_evd_m(Complex* A, int N, Complex* u, Complex* s);
void jacobi_evd_n(Complex* A, int N, Complex* u, Complex* s);
