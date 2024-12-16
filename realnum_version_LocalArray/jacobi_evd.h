#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double* generate_matrix(int n);
void print_matrix(double*A, int N);
void make_identity(double* u, int N);
void transpose(double* A, double* A_t, int N);
void matmul(double*A, double* B, double* C, int N);
void jacobi_rotate(double *A, double *U, int N, int k, int l);
void find_maxDiagOff(double* A, int N, int* p, int *q, double * value);
double cal_sumDiagOff(double *A, int N);
void jacobi_evd(double* A, int N, double* u, double* s);
void jacobi_evd_m(double* A, int N, double* u, double* s);
void jacobi_evd_n(double* A, int N, double* u, double* s);
void jacobi_test(double* orgA, int N, double* u, double* s);
