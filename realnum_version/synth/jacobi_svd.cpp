#include "jacobi_svd.h"
#include <hls_stream.h>
#include "ap_axi_sdata.h"
#include <hls_math.h>
#include <ap_fixed.h>

#define M 10
#define N 20

typedef hls::axis<double, 0, 0, 1> axis_double;

double A_t[N*M];
double AA_t[M*M]; 
double A_tA[N*N]; 
double s_m[M*M];
double s_n[N*N];
double u_t[M*M];
double v_t[N*N];
double us[M*N];
double newU[M*M];
double tmp[M*N]; 
double uSv_t[M*N]; 
double vs_t[N*M];
double newV[N*N];
double newV_t[N*N];
double s_t[N*M];
        
void jacobi_svd(hls::stream<axis_double>& A_in_stream, hls::stream<axis_double>& U_out_stream, hls::stream<axis_double>& V_out_stream, hls::stream<axis_double>& S_out_stream, int m, int n) {

	#pragma HLS INTERFACE axis port=A_in_stream
	#pragma HLS INTERFACE axis port=U_out_stream
	#pragma HLS INTERFACE axis port=V_out_stream
	#pragma HLS INTERFACE axis port=S_out_stream
	#pragma HLS INTERFACE s_axilite port=m bundle=control
	#pragma HLS INTERFACE s_axilite port=n bundle=control
	#pragma HLS INTERFACE s_axilite port=return bundle=control

	// A : m by n
    // u: m by m 
    // v: n by n 
    // s: m by n 

    double A[M * N];
    double U[M * M];
    double V[N * N];
    double S[M * N];

    // Read input matrix from AXI stream
    for (int i = 0; i < m * n; i++) {
        axis_double tmp = A_in_stream.read();
        A[i] = tmp.data;
    }

    /* ----- STEP 0 -----*/
    transpose_mxn(A, A_t, m, n);
    matmul_mxn(A, A_t, AA_t, m, n, n, m);
    matmul_mxn(A_t, A, A_tA, n, m, m, n);
    
    /* ----- STEP 1 -----*/
    // Initialize u
    make_identity(U, m);

    // Initialize s_m : (m x m diagonal matrix)
    make_identity(s_m, m);

    // RUN Jacobi EVD for AA_t : AA_t = (U)*(s_m)*(U^t)
    jacobi_evd_m(AA_t, m, U, s_m);
    eigen_sort(U, s_m, m); // Sort decreasing diagonal values of s_m & corresponding eigen vector u_i

    /* ----- STEP 2 -----*/
    // Initialize v
    make_identity(V, n);
    
    // Initialize s_n : (n x n diagonal matrix)
    make_identity(s_n, n); 
    
    // RUN Jacobi EVD for A_tA : A_tA = (V)*(s_n)*(V^t)
    jacobi_evd_n(A_tA, n, V, s_n);
    eigen_sort(V, s_n, n); // Sort decreasing diagonal values of s_n & corresponding eigen vector v_i

    /* ----- STEP 3 -----*/
    // Update Scale Matrix called SIGMA (m x n matrix s)
    make_scale_mxn(s_m, s_n, S, m, n);
    printf("------Scale matrix SIGMA-----\n");
    print_mxn_matrix(S, m, n);
    
    /*----- STEP 4 -----*/
    // Jacobi SVD  

    transpose(U, u_t, m);
    transpose(V, v_t, n);
    
    if (m <= n) {
        // (m < n) Get U without calculating eig(AA^T)
        matmul_mxn(A, V, us, m, n, n, n);
        get_u_from_us(us, S, newU, m, n);
        matmul_mxn(newU, S, tmp, m, m, m, n);
        matmul_mxn(tmp, v_t, uSv_t, m, n, n, n);

        // Return U
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                U[i * m + j] = newU[i * m + j];
            }
        }
    } else {
        // m > n, Get V without calculating eig(A^T A)
        matmul_mxn(A_t, U, vs_t, n, m, m, m);
        transpose_mxn(S, s_t, m, n);
              
        get_v_from_vs_t(vs_t, s_t, newV, n, m); 
        
        transpose_mxn(newV, newV_t, n, n);
        matmul_mxn(U, S, tmp, m, m, m, n);
        matmul_mxn(tmp, newV_t, uSv_t, m, n, n, n);

        // Return V
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                V[i * n + j] = newV[i * n + j];
            }
        }
        // U remains the same
    }

    // Write U, S, V matrices to output AXI stream
    for (int i = 0; i < m * m; i++) {
        axis_double tmp;
        tmp.data = U[i];
        U_out_stream.write(tmp);
    }

    for (int i = 0; i < n * n; i++) {
        axis_double tmp;
        tmp.data = V[i];
        V_out_stream.write(tmp);
    }

    for (int i = 0; i < m * n; i++) {
        axis_double tmp;
        tmp.data = S[i];
        S_out_stream.write(tmp);
    }
}
