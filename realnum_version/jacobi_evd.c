#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <math.h> 
#include "jacobi_evd.h"

#define m 10
#define n 20


void print_matrix(double*A, int N){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            printf("%f " , A[i*N+j]); 
        }
        printf("\n");
    }
    printf("\n");
}

void make_identity(double* u, int N){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            if(i==j){
                u[i*N+j]= 1; 
            }
            else{
                u[i*N+j] = 0;
            }
        }
    }
}


void transpose(double* A, double* A_t, int N){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            A_t[j*N+i] = A[i*N+j];
        }
    }
 
}

void matmul(double*A, double* B, double* C, int N){
    // Assume A , B, C : N by N square matrix
    
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            C[i*N +j] = 0; 
            for(int k=0; k<N; k++){
                C[i*N + j] += A[i*N+k] * B[k*N+j]; 
            }
        }
    }
    
}
void jacobi_rotate(double *A, double *U, int N, int k, int l){
    double t, c, s; 

    double a_kk = A[k*N + k];
    double a_ll = A[l*N + l];
    double a_kl = A[k*N + l];

    double theta = 0.5* atan2(2*a_kl, a_ll-a_kk); 

    t = tan(theta); 
    c = 1/(sqrt(1+t*t)); 
    s = c*t; 


    //update matrix A
    
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

void find_maxDiagOff(double* A, int N, int* p, int *q, double * value){

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

double cal_sumDiagOff(double *A, int N){
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

void jacobi_evd_m(double* A, int N, double* u, double* s){
    
    //double *A_t, *AA_t; 
    //A_t = (double*)malloc(sizeof(double)*N*N); 

    static double A_t[m*m];      

    //AA_t = (double*)malloc(sizeof(double)*N*N);
    static double AA_t[m*m];

    transpose(A, A_t, N);
    matmul(A, A_t, AA_t, N); 

    //double* diagoffSum = (double*)malloc(sizeof(double)); 
    double diagoffSum[1];

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

void jacobi_evd(double* A, int N, double* u, double* s){
    
    //double *A_t, *AA_t; 
    //A_t = (double*)malloc(sizeof(double)*N*N); 

    double A_t[N*N];      

    //AA_t = (double*)malloc(sizeof(double)*N*N);
    double AA_t[N*N];

    transpose(A, A_t, N);
    matmul(A, A_t, AA_t, N); 

    //double* diagoffSum = (double*)malloc(sizeof(double)); 
    double diagoffSum[1];

    int rot_num = 0 ; 
    while(cal_sumDiagOff(A, N)>1e-10){
        
        //int* p = (int*)malloc(sizeof(int)); 
        int p_ = -1;
        int* p = &p_; 

        //int *q = (int*)malloc(sizeof(int)); ; 
        int q_ = 0;
        int * q = &q_; 

        //double *value = (double*)malloc(sizeof(double)); ; 
        double* value;
    
        find_maxDiagOff(A, N, p, q, value); 

        jacobi_rotate(A, u,N,*p,*q); 
        
        rot_num ++; 
    }
    
    //update s 
    for(int i=0; i<N; i++){
        s[i*N+i] = A[i*N+i]; 
    }


}

