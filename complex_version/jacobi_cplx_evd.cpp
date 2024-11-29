#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <math.h> 
#include "jacobi_cplx_evd.h"
#include <iostream> 

double complex_norm(Complex z){
    return sqrt((z.real() * z.real()) + (z.imag()*z.imag())); 
}

void print_matrix(Complex*A, int N){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            //printf("%f " , A[i*N+j]); 
            std::cout<< A[i*N+j] << " ";
        }
        printf("\n");
    }
    printf("\n");
}

void make_identity(Complex* u, int N){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            if(i==j){
                //u[i*N+j]= 1; 
                u[i*N+j]= Complex(1,0); 
            }
            else{
                //u[i*N+j] = 0;
                u[i*N+j]= Complex(0,0); 
            }
        }
    }
}


void transpose(Complex* A, Complex* A_t, int N){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            A_t[j*N+i] = A[i*N+j];
        }
    }
 
}

void matmul(Complex* A, Complex*  B, Complex* C, int N){
    // Assume A , B, C : N by N square matrix
    
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            //C[i*N +j] = 0; 
            C[i*N +j] = Complex(0,0); 
            for(int k=0; k<N; k++){
                C[i*N + j] += A[i*N+k] * B[k*N+j]; 
            }
        }
    }
    
}
// Should Change to get theta in Complex version 
void jacobi_rotate(Complex *A, Complex *U, int N, int k, int l){
    double t, c, s; 

    Complex a_kk = A[k*N + k];
    Complex a_ll = A[l*N + l];
    Complex a_kl = A[k*N + l];

    // Should Change to get theta in Complex version 
    double theta = 0.5* atan2(2*a_kl, a_ll-a_kk); 

    t = tan(theta); 
    c = 1/(sqrt(1+t*t)); 
    s = c*t; 


    //update matrix A
    // Should Change to get theta in Complex version 
    
    for(int h=0; h<N; h++){
        double a_hk = A[N*h+k];
        double a_hl = A[N*h+l];

        A[N*h+k] = A[N*k+h] = c*a_hk - s*a_hl;
        A[N*h+l] = A[N*l+h] = c*a_hl + s*a_hk; 

    }
    A[N*k+l] = A[N*l+k] = 0; 
    
    A[N*k+k] = c*c*a_kk + s*s*a_ll - 2*s*c*a_kl;
    A[N*l+l] = s*s*a_kk + c*c*a_ll + 2*s*c*a_kl; 
    
    //update U (sample)
    double u_kk = U[k*N+k];
    double u_kl = U[k*N+l];
    double u_lk = U[l*N+k];
    double u_ll = U[l*N+l]; 

    for(int i=0; i<N; i++){

        double u_ik = U[i*N + k];
        double u_il = U[i*N + l];

        U[i*N+k]= c*u_ik - s*u_il; 
        U[i*N+l]= s*u_ik + c*u_il;
    }


}
// Should Change to get theta in Complex version 
void find_maxDiagOff(Complex* A, int N, int* p, int *q, double* value){

    *value = -1;
    *p=-1; 
    *q=-1; 

    for (int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            
            double tmp = A[i*N+j]; 
            
            if(tmp<0){
                tmp = (-1)*tmp;
            }

            if(*value < tmp){
                *value = tmp;
                *p = i;  
                *q = j; 
            }
            
        }
    }
}

// Should Change to get theta in Complex version 
double cal_sumDiagOff(Complex *A, int N){
    double ret=0; 
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            if (A[i*N+j]<0.0){
                ret += (-1)*A[i*N+j];
            }
            else{
                ret += A[i*N+j];
            }
        }
    }
    return ret; 
}

// Should Change to get theta in Complex version 
void jacobi_evd_m(Complex* A, int N, Complex* u, Complex* s){
    
    //double *A_t, *AA_t; 
    //A_t = (double*)malloc(sizeof(double)*N*N); 

    static Complex A_t[m*m];      

    //AA_t = (double*)malloc(sizeof(double)*N*N);
    static Complex AA_t[m*m];

    transpose(A, A_t, N);
    matmul(A, A_t, AA_t, N); 

    //double* diagoffSum = (double*)malloc(sizeof(double)); 
    //double diagoffSum[1];

    int rot_num = 0 ; 
    while(cal_sumDiagOff(A, N)>1e-10){
        
        //int* p = (int*)malloc(sizeof(int)); 
        //int p_ = -1;
        //int* p = &p_; 
        int p = -1 ; 

        //int *q = (int*)malloc(sizeof(int)); ; 
        //int q_ = 0;
        //int * q = &q_; 
        int q = -1; 

        //double *value = (double*)malloc(sizeof(double)); ; 
        //double* value;
        double value= -1; 
        find_maxDiagOff(A, N, &p, &q, &value); 

        jacobi_rotate(A, u,N,p,q); 
        
        rot_num ++; 
    }
    
    //update s 
    for(int i=0; i<N; i++){
        s[i*N+i] = A[i*N+i]; 
    }


}


void jacobi_evd_n(double* A, int N, double* u, double* s){
    
    //double *A_t, *AA_t; 
    //A_t = (double*)malloc(sizeof(double)*N*N); 

    static double A_t[n*n];      

    //AA_t = (double*)malloc(sizeof(double)*N*N);
    static double AA_t[n*n];

    transpose(A, A_t, N);
    matmul(A, A_t, AA_t, N); 

    //double* diagoffSum = (double*)malloc(sizeof(double)); 
    double diagoffSum[1];

    int rot_num = 0 ; 
    while(cal_sumDiagOff(A, N)>1e-10){
        
        //int* p = (int*)malloc(sizeof(int)); 
        //int p_ = -1;
        //int* p = &p_; 
        int p = -1; 

        //int *q = (int*)malloc(sizeof(int)); ; 
        //int q_ = 0;
        //int * q = &q_; 
        int q = -1; 

        //double *value = (double*)malloc(sizeof(double)); ; 
        //double* value;
        double value = 0;
    
        find_maxDiagOff(A, N, &p, &q, &value); 

        jacobi_rotate(A, u,N,p,q); 
        
        rot_num ++; 
    }
    
    //update s 
    for(int i=0; i<N; i++){
        s[i*N+i] = A[i*N+i]; 
    }


}
