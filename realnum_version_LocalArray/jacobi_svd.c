//#include "mat_utils.h"
//#include "jacobi_evd.h"
#include "jacobi_svd.h"

#define M 10
#define N 20

double A_t[N*M];
double AA_t[M*M]; 
double A_tA[N*N]; 
double s_m[M*M];
double s_n[N*N] ; 
double u_t[M*M] ;
double v_t[N*N] ;
double us[M*N];
double newU[M*M];
double tmp[M*N]; 
double uSv_t[M*N]; 
double vs_t[N*M] ; 
double newV[N*N] ; 
double newV_t[N*N] ;
double s_t[N*M] ; 
//double tmp[m*n]; 
//double uSv_t[m*n]; 
        

void jacobi_svd(double* A, double* u, double*v, double*s, int m, int n){

    // A : m by n 
    // u: m by m 
    // v: n by n 
    // s: m by n 

    // jacobi_svd returns u, v, s such as A = (u)(s)(v^t)
    
    /* ----- STEP 0 -----*/
    
    // get A [m by n] matrix 
    //printf("-----A-----\n");
    //print_mxn_matrix(A, m, n);

    // get A_t, AA_t, A_tA matrix 

    // double A_t[n*m];
    // double AA_t[m*m]; 
    // double A_tA[n*n]; 
    transpose_mxn(A, A_t, m,n);     
    matmul_mxn(A, A_t, AA_t, m,n,n,m); 
    matmul_mxn(A_t, A, A_tA, n,m,m,n);
    
    /* ----- STEP 1 -----*/

    // copy AA_t  
    
    //double orgAA_t[m*m] ; 
    //memcpy(orgAA_t, AA_t, sizeof(double)*m*m); 

    // initialize u  
    make_identity(u, m); 

    // initialize s_m : (mxm diagonal matrix) 
    //double s_m[m*m];
    make_identity(s_m, m);

    // RUN jacobi EVD for AA_t : AA_t = (u)*(s_m)*(u^t)
    jacobi_evd_m(AA_t, m, u, s_m);
    //jacobi_evd_m(AA_t, m, u, s_m);
    eigen_sort(u, s_m, m); // Sort decreasingly diagonal values of s_m & corrensponded eigen vector u_i   
    //printf("-----------EVD test : U--------\n");
    //print_mxn_matrix(u, m,m);

    /* ----- STEP 2 -----*/
    
    // copy A_tA 
    //double orgA_tA[n*n]; 
    //memcpy(orgA_tA, A_tA, sizeof(double)*n*n); 

    // initialize v 
    make_identity(v, n); 
    
    //initialize s_n : (nxn diagonal matrix)
    //double s_n[n*n] ; 
    make_identity(s_n, n); 
    
    // RUN jacobi EVD for A_tA : A_tA = (v)*(s_n)*(v^t) 
    jacobi_evd_n(A_tA, n, v, s_n); 
    eigen_sort(v, s_n, n); // Sort decreasingly diagonal values of s_n & corrensponded eigen vector v_i   
    //printf("-----------EVD test : V--------\n");
    //print_mxn_matrix(v, n,n);

    /* ----- STEP 3 -----*/
    // update Scale Matrix called SIGMA mxn matrix s 
    
    // double s_mxn[m*n];
    make_scale_mxn(s_m, s_n, s, m, n) ;
    //printf("------scale matrix SIGMA-----\n");
    //print_mxn_matrix(s, m,n);
    
    /*----- STEP 4 -----*/
    // Jacobi SVD  
    
    // get u_t and v_t 
    // double u_t[m*m] ;
    // double v_t[n*n] ;

    transpose(u, u_t, m);     
    transpose(v, v_t, n); 
    
    if (m<=n){
        // (m<n) get U without calculating eig(AAt); 
        // double us[m*n];
        // double newU[m*m];
        
        matmul_mxn(A, v, us, m, n, n, n);
        /*
        printf("----A----\n");
        print_mxn_matrix(A, m,n);
        
        printf("----v----\n");
        print_mxn_matrix(v, n,n);
        printf("----us----\n");
        print_mxn_matrix(us, m,n);
        printf("----s1----\n");
        print_mxn_matrix(s, m,n);
        */
        
        get_u_from_us(us, s, newU, m, n); 
        //printf("----s2----\n");
        //print_mxn_matrix(s, m,n);
        

        
        //printf("----newU----\n");
        //print_mxn_matrix(newU, m,m);

        //printf("----V^t----\n");
        //print_mxn_matrix(v_t, m,m);

        // double tmp[m*n]; 
        // double uSv_t[m*n]; 
        matmul_mxn(newU, s, tmp, m, m, m, n);
        matmul_mxn(tmp, v_t, uSv_t,m,n,n,n); 
        //printf("---Jacobi SVD Result: uSv_t---\n"); 
        //print_mxn_matrix(uSv_t, m, n);
        
        //return u, v 
        //memcpy(u, newU, sizeof(double)*m*m); 
        // replace memcpy 

        for(int i=0; i<m; i++){
            for(int j=0; j<m; j++){
                u[i*m+j] = newU[i*m+j];
            }
        }


    }
    else{ // m> n
        // get V without calculating eig(AtA); 
   
        // double vs_t[n*m] ; 
        // double newV[n*n] ; 
        // double newV_t[n*n] ;
        // double s_t[n*m] ; 
        
        matmul_mxn(A_t, u, vs_t, n, m, m, m); 
        transpose_mxn(s, s_t, m, n);
              
        get_v_from_vs_t(vs_t, s_t, newV, n, m); 
        
        transpose_mxn(newV, newV_t, n,n); 

        //printf("----U----\n");
        //print_mxn_matrix(u, m,m);

        //printf("----V^t----\n");
        //print_mxn_matrix(newV_t, n,n); 
        
        //double tmp[m*n]; 
        matmul_mxn(u, s, tmp, m, m, m, n);
        //double uSv_t[m*n]; 
        matmul_mxn(tmp, newV_t, uSv_t,m,n,n,n); 
        //printf("---Jacobi SVD Result: uSv_t---\n"); 
        //print_mxn_matrix(uSv_t, m, n);

        //return u, v 
        //memcpy(v, newV, sizeof(double)*n*n);

        //replace memcpy 
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                v[i*n+j] = newV[i*n+j]; 
            }
        }

        // u is same 
    }

}
