#include "mat_utils_cplx.h"


void print_mxn_matrix(Complex* A, int m, int n){

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            //printf("%f ",A[i*n+j]); 
            std::cout << A[i*n+j] << " "; 
        }
        std::cout << std::endl; 
    }
    std::cout << std::endl; 

}

void matmul_mxn(Complex* A, Complex* B, Complex*C, int rowA, int colA, int rowB, int colB){

    // A : rowA X colA 
    // B : rowB X colB (colA == rowB)
    // C (= AXB) : rowA X colB 

    if(colA != rowB){
        //printf("col# of A != row#B\n"); 
        std::cout << "col# of A != row#B" << std::endl;
        return; 
    }

    for(int i=0; i<rowA; i++){
        for(int j=0; j<colB; j++){
            Complex AikBkj = Complex(0,0) ;
            for(int k=0; k<colA; k++){
                AikBkj += A[i*colA + k] * B[k*colB+ j];
                //printf("i: %d j: %d : %d * %d\n",i,j, A[i*rowA+k], A[k*rowB+j] );
            }
            C[i*colB + j] = AikBkj;
            //printf("C_{%d}_{%d}: %f\n", i, j, AikBkj);
        }
    }
} 

void transpose_mxn(Complex* A, Complex* A_t, int m, int n){
    // A: m x n matrix 
    // A_t : n x m matrix 

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            A_t[j*m+i] = A[i*n+j]; 
        }
    }
}

void hermitian_mxn(Complex*A, Complex* A_H, int m, int n){
    // Conjugate Transpose 
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            A_H[j*m+i] = conj(A[i*n+j]);
        }
    }    
}

void make_scale_mxn(Complex* s_m, Complex* s_n, Complex* s_mxn, int m, int n){
    
    //int s = -1; 
    if (m>n){
        //s = n;

        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                if(i==j){
                    
                    //s_mxn[i*n +j] = sqrt(fabs(s_n[i*n+j])); 
                    s_mxn[i*n +j] = Complex(sqrt(complex_norm(s_n[i*n+j])), 0); 
                }
                else{
                    //s_mxn[i*n + j] = 0; 
                    s_mxn[i*n +j] = Complex(0,0);
                }
            }
        } 
    }
    else{ // m<=n
        //s = m; 

        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                if(i==j){
                    s_mxn[i*n+j] = Complex(sqrt(complex_norm(s_m[i*m+j])), 0); 
                    //s_mxn[i*n+j] = sqrt(fabs(s_m[i*m+j])); 
                }
                else{
                    s_mxn[i*n+j]= Complex(0,0); 
                }
            }
        }
    }

} 

void swap_vectors(int mth, int nth, Complex* EV, int N){

    for(int i=0; i<N; i++){
        // EV[i][col_i] <--> EV[i][col_j] 
        Complex ev_im = EV[i*N + mth];
        Complex ev_in = EV[i*N + nth];

        EV[i*N + mth] = ev_in; 
        EV[i*N + nth] = ev_im; 
    }
}

void eigen_sort(Complex* EigenVectors, Complex* s, int m){
    
    // bubble sort (descending order)

    for(int i=0; i<m-1; i++){

        for(int j=0; j<m-i-1; j++){

            Complex sjj = s[j*m + j]; // Sjj
            Complex sjj_next = s[(j+1)*m + j+1]; //Sj+1, j+1

            if(complex_norm(sjj) < complex_norm(sjj_next)){
                
                s[(j+1)*m + (j+1)] = sjj;  
                s[j*m + j] = sjj_next;

                swap_vectors(j, j+1, EigenVectors, m);  
            }
        }
    }

}

void get_u_from_us(Complex* us, Complex* s_mxn, Complex* newU, int m, int n){

    //us : mxn  (us==Av)  
    // s_mxn : mxn 
    // newU : mxm

    if(m<=n){
        for(int j=0; j<m; j++){ // jth column 
            for(int i=0; i<m; i++){ // ith row 
                if (complex_norm(s_mxn[j*n+j])< 10e-7){
                    s_mxn[j*n+j] = 10e-8; 
                }
                newU[i*m + j] = us[i*n+j] / s_mxn[j*n+j];    
            }
        }
    }    

}

void get_v_from_vs_t(Complex* vs_t, Complex* s_mxn_t, Complex* newV, int n, int m){
    
    // (v_t)(s_t) : nxm (uA == Sv --> (Sv)^t = (uA)^t) 
    // s_mxn_t : nxm  
    // newV_t : nxn 

    if(m>n){
        for(int j=0; j<n; j++){ //j th column 
            for(int i=0; i<n; i++){ // ith row 

                if (complex_norm(s_mxn_t[j*m+j])< 10e-7){
                    s_mxn_t[j*m+j] = 10e-8; 
                }
                newV[n*i + j] = vs_t[i*m+j] / s_mxn_t[j*m+j] ; 
            }
        }
    }

}