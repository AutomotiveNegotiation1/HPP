#pragma once 

#include <stdlib.h> 
#include <stdio.h> 
#include <string.h>
#include <math.h>
#include "jacobi_config.h"
#include "jacobi_type.h"
#include <iostream>
#include "jacobi_cplx_evd.h"

void print_mxn_matrix(Complex* A, int rowA, int colA);
void matmul_mxn(Complex* A, Complex* B, Complex*C, int rowA, int colA, int rowB, int colB);
void transpose_mxn(Complex* A, Complex* A_t, int rowA, int colA);
void hermitian_mxn(Complex*A, Complex* A_H, int rowA, int colA);

void make_scale_mxn(Complex* s_m, Complex* s_n, Complex* s_mxn, int rowS, int colS);
void eigen_sort(Complex* EigenVectors, Complex* s, int N);
void swap_vectors(int mth, int nth, Complex* EV, int N);
void get_u_from_us(Complex* us, Complex* s_mxn, Complex* newU, int rowS, int colS);
void get_v_from_vs_t(Complex* v_ts_t, Complex* s_mxn_t, Complex* newV_t, int rowSt, int colSt);