#pragma once 

#include "mat_utils.h"
#include "jacobi_evd.h"


void jacobi_svd(double* A, double* u, double*v, double*s, int m, int n);