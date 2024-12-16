#include <stdlib.h> 
#include <stdio.h> 
#include <string.h>
#include <math.h>


void print_mxn_matrix(double*A, int m, int n);
void matmul_mxn(double* A, double* B, double*C, int rowA, int colA, int rowB, int colB);
void transpose_mxn(double* A, double* A_t, int m, int n);

void make_scale_mxn(double* s_m, double* s_n, double* s_mxn, int m, int n);
void eigen_sort(double* EigenVectors, double* s, int m);
void swap_vectors(int mth, int nth, double* EV, int N);
void get_u_from_us(double* us, double* s_mxn, double* newU, int m, int n);
void get_v_from_vs_t(double* v_ts_t, double* s_mxn_t, double* newV_t, int n, int m);