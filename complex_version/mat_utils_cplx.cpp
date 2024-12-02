#include "mat_utils_cplx.h"


void print_mxn_matrix(Complex* A, int rowA, int colA){

    for(int i=0; i<rowA; i++){
        for(int j=0; j<colA; j++){
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

void transpose_mxn(Complex* A, Complex* A_t, int rowA, int colA){
    // A: m x n matrix 
    // A_t : n x m matrix 

    for(int i=0; i<rowA; i++){
        for(int j=0; j<colA; j++){
            A_t[j*rowA+i] = A[i*colA+j]; 
        }
    }
}

void hermitian_mxn(Complex*A, Complex* A_H, int rowA, int colA){
    // Conjugate Transpose 
    for(int i=0; i<rowA; i++){
        for(int j=0; j<colA; j++){
            A_H[j*rowA+i] = conj(A[i*colA+j]);
        }
    }    
}

void make_scale_mxn(Complex* s_m, Complex* s_n, Complex* s_mxn, int rowS, int colS){
    
    //int s = -1; 
    
    // m->rowS 
    // n->colS 

    if (rowS>colS){
        //s = n;

        for(int i=0; i<rowS; i++){
            for(int j=0; j<colS; j++){
                if(i==j){
                    
                    //s_mxn[i*n +j] = sqrt(fabs(s_n[i*n+j])); 
                    s_mxn[i*colS +j] = Complex(sqrt(complex_norm(s_n[i*colS+j])), 0); 
                }
                else{
                    //s_mxn[i*n + j] = 0; 
                    s_mxn[i*colS +j] = Complex(0,0);
                }
            }
        } 
    }
    else{ // m<=n
        //s = m; 

        for(int i=0; i<rowS; i++){
            for(int j=0; j<colS; j++){
                if(i==j){
                    s_mxn[i*colS+j] = Complex(sqrt(complex_norm(s_m[i*rowS+j])), 0); 
                    //s_mxn[i*n+j] = sqrt(fabs(s_m[i*m+j])); 
                }
                else{
                    s_mxn[i*colS+j]= Complex(0,0); 
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

void eigen_sort(Complex* EigenVectors, Complex* s, int N){
    
    // bubble sort (descending order)

    // m -> N 
    for(int i=0; i<N-1; i++){

        for(int j=0; j<N-i-1; j++){

            Complex sjj = s[j*N + j]; // Sjj
            Complex sjj_next = s[(j+1)*N + j+1]; //Sj+1, j+1

            if(complex_norm(sjj) < complex_norm(sjj_next)){
                
                s[(j+1)*N + (j+1)] = sjj;  
                s[j*N + j] = sjj_next;

                swap_vectors(j, j+1, EigenVectors, N);  
            }
        }
    }

}

void get_u_from_us(Complex* us, Complex* s_mxn, Complex* newU, int rowS, int colS){

    //us : mxn  (us==Av)  
    // s_mxn : mxn 
    // newU : mxm

    //m -> rowS 
    //n -> colS 

    if(rowS<=colS){
        for(int j=0; j<rowS; j++){ // jth column 
            for(int i=0; i<rowS; i++){ // ith row 
                if (complex_norm(s_mxn[j*colS+j])< 10e-7){
                    s_mxn[j*colS+j] = 10e-8; 
                }
                newU[i*rowS + j] = us[i*colS+j] / s_mxn[j*colS+j];    
            }
        }
    }    

}

void get_v_from_vs_t(Complex* vs_t, Complex* s_mxn_t, Complex* newV, int rowSt, int colSt){
    
    // (v_t)(s_t) : nxm (uA == Sv --> (Sv)^t = (uA)^t) 
    // s_mxn_t : nxm  
    // newV_t : nxn 

    // m -> colSt
    // n --> rowSt

    if(colSt>rowSt){
        for(int j=0; j<rowSt; j++){ //j th column 
            for(int i=0; i<rowSt; i++){ // ith row 

                if (complex_norm(s_mxn_t[j*colSt+j])< 10e-7){
                    s_mxn_t[j*colSt+j] = 10e-8; 
                }
                newV[rowSt*i + j] = vs_t[i*colSt+j] / s_mxn_t[j*colSt+j] ; 
            }
        }
    }

}