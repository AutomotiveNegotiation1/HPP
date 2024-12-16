#include "mat_utils.h"


#define M_max 10
#define N_max 20


void print_mxn_matrix(double*A, int m, int n){

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            printf("%f ",A[i*n+j]); 
        }
        printf("\n");
    }
    printf("\n");

}

void matmul_mxn(double* A, double* B, double*C, int rowA, int colA, int rowB, int colB){

    // A : rowA X colA 
    // B : rowB X colB (colA == rowB)
    // C (= AXB) : rowA X colB 

    if(colA != rowB){
        printf("col# of A != row#B\n"); 
        return; 
    }

	#pragma HLS PIPELINE II=1
    for(int i=0; i<rowA; i++){
        for(int j=0; j<colB; j++){
            double AikBkj = 0 ;
            for(int k=0; k<colA; k++){
				#pragma HLS UNROLL factor=4
                AikBkj += A[i*colA + k] * B[k*colB+ j];
                //printf("i: %d j: %d : %d * %d\n",i,j, A[i*rowA+k], A[k*rowB+j] );
            }
            C[i*colB + j] = AikBkj;
            //printf("C_{%d}_{%d}: %f\n", i, j, AikBkj);
        }
    }
} 

void transpose_mxn(double* A, double* A_t, int m, int n){
    // A: m x n matrix 
    // A_t : n x m matrix 

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            A_t[j*m+i] = A[i*n+j]; 
        }
    }
}

void make_scale_mxn(double* s_m, double* s_n, double* s_mxn, int m, int n){
    
    //int s = -1; 
    if (m>n){
        //s = n;

        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                if(i==j){
                    
                    s_mxn[i*n +j] = sqrt(fabs(s_n[i*n+j])); 
                }
                else{
                    s_mxn[i*n + j] = 0; 
                }
            }
        } 
    }
    else{ // m<=n
        //s = m; 

        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                if(i==j){
                    
                    s_mxn[i*n+j] = sqrt(fabs(s_m[i*m+j])); 
                }
                else{
                    s_mxn[i*n+j]= 0; 
                }
            }
        }
    }

} 

void swap_vectors(int mth, int nth, double* EV, int N){
	#pragma HLS PIPELINE
    for(int i=0; i<N; i++){
		#pragma HLS UNROLL factor=2
        // EV[i][col_i] <--> EV[i][col_j] 
//        double ev_im = EV[i*N + mth];
//        double ev_in = EV[i*N + nth];
//        EV[i*N + mth] = ev_in;
//        EV[i*N + nth] = ev_im;
        double temp = EV[i * N + mth];
        EV[i * N + mth] = EV[i * N + nth];
        EV[i * N + nth] = temp;
    }
}

void eigen_sort(double*EigenVectors, double* s, int m){
    
    // bubble sort (descending order)
    for(int i=0; i<m-1; i++){
		#pragma HLS PIPELINE
        for(int j=0; j<m-i-1; j++){
			#pragma HLS UNROLL factor=2
            double sjj = s[j*m + j]; // Sjj
            double sjj_next = s[(j+1)*m + j+1]; //Sj+1, j+1

            if(sjj < sjj_next){
                
//                s[(j+1)*m + (j+1)] = sjj;
//                s[j*m + j] = sjj_next;
                double temp = sjj;
                s[j * m + j] = sjj_next;
                s[(j + 1) * m + (j + 1)] = temp;

                swap_vectors(j, j+1, EigenVectors, m);  
            }
        }
    }
}

//void get_u_from_us(double* us, double* s_mxn, double* newU, int m, int n){
//
//    //us : mxn  (us==Av)
//    // s_mxn : mxn
//    // newU : mxm
//
//    if(m<=n){
//        for(int j=0; j<m; j++){ // jth column
//            for(int i=0; i<m; i++){ // ith row
//                if (fabs(s_mxn[j*n+j])< 10e-7){
//                    s_mxn[j*n+j] = 10e-8;
//                }
//                newU[i*m + j] = us[i*n+j] / s_mxn[j*n+j];
//            }
//        }
//    }
//
//}

void get_u_from_us(double* us, double* s_mxn, double* newU, int m, int n) {
    #pragma HLS ARRAY_PARTITION variable=us cyclic factor=4 dim=1
    #pragma HLS ARRAY_PARTITION variable=s_mxn cyclic factor=4 dim=1
    #pragma HLS ARRAY_PARTITION variable=newU cyclic factor=4 dim=1

    if (m <= n) {
        double s_values[M_max];
		#pragma HLS ARRAY_PARTITION variable=s_values complete dim=1

        for (int j = 0; j < m; j++) {
            double value = fabs(s_mxn[j * n + j]);
            s_values[j] = (value < 10e-7) ? 10e-8 : value;
        }

        for (int j = 0; j < m; j++) { // jth column
            #pragma HLS UNROLL factor=4
            for (int i = 0; i < m; i++) { // ith row
                newU[i * m + j] = us[i * n + j] / s_values[j];
            }
        }
    }
}

//void get_v_from_vs_t(double* vs_t, double* s_mxn_t, double* newV, int n, int m){
//
//    // (v_t)(s_t) : nxm (uA == Sv --> (Sv)^t = (uA)^t)
//    // s_mxn_t : nxm
//    // newV_t : nxn
//
//    if(m>n){
//        for(int j=0; j<n; j++){ //j th column
//            for(int i=0; i<n; i++){ // ith row
//
//                if (fabs(s_mxn_t[j*m+j])< 10e-7){
//                    s_mxn_t[j*m+j] = 10e-8;
//                }
//                newV[n*i + j] = vs_t[i*m+j] / s_mxn_t[j*m+j] ;
//            }
//        }
//    }
//
//}


void get_v_from_vs_t(double* vs_t, double* s_mxn_t, double* newV, int n, int m) {
    #pragma HLS ARRAY_PARTITION variable=vs_t cyclic factor=4 dim=1
    #pragma HLS ARRAY_PARTITION variable=s_mxn_t cyclic factor=4 dim=1
    #pragma HLS ARRAY_PARTITION variable=newV cyclic factor=4 dim=1

    if (m > n) {
        double s_values[N_max];
		#pragma HLS ARRAY_PARTITION variable=s_values complete dim=1

        for (int j = 0; j < n; j++) {
            double value = fabs(s_mxn_t[j * m + j]);
            s_values[j] = (value < 10e-7) ? 10e-8 : value;
        }


        for (int j = 0; j < n; j++) { // j th column
            #pragma HLS UNROLL factor=4
            for (int i = 0; i < n; i++) { // i th row
                newV[n * i + j] = vs_t[i * m + j] / s_values[j];
            }
        }
    }
}
