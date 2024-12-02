//#include "mat_utils.h"
#include "jacobi_cplx_evd.h"
#include "jacobi_config.h"
#include "jacobi_cplx_svd.h"

        

void jacobi_svd(Complex* A, Complex* u, Complex* v, Complex* s){

    // A : m by n 
    // u: m by m 
    // v: n by n 
    // s: m by n 

    // jacobi_svd returns u, v, s such as A = (u)(s)(v^t)


    Complex A_H[n*m];
    Complex AA_H[m*m]; 
    Complex A_HA[n*n];

    Complex s_m[m*m];
    Complex s_n[n*n] ; 

    Complex u_H[m*m] ;
    Complex v_H[n*n] ;

    Complex us[m*n];
    Complex newU[m*m];
    Complex tmp[m*n]; 
    Complex uSv_H[m*n]; 
    Complex vs_H[n*m] ; 

    Complex newV[n*n] ; 
    Complex newV_H[n*n] ;
    Complex s_H[n*m] ; 
    
    /* ----- STEP 0 -----*/
    
    // get A [m by n] matrix 
    std::cout << "------A-----"<< std::endl;
    print_mxn_matrix(A, m, n);

    // get A_H, AA_H, A_HA matrix 
    //transpose_mxn(A, A_t, m,n); 
    hermitian_mxn(A, A_H, m, n);

    matmul_mxn(A, A_H, AA_H, m,n,n,m); 
    matmul_mxn(A_H, A, A_HA, n,m,m,n);
    
    /* ----- STEP 1 -----*/

    // initialize u  
    make_identity(u, m); 

    // initialize s_m : (mxm diagonal matrix) 
    //double s_m[m*m];
    make_identity(s_m, m);

    // RUN jacobi EVD for AA_t : AA_t = (u)*(s_m)*(u^t)
    jacobi_evd_m(AA_H, m, u, s_m);
    
    eigen_sort(u, s_m, m); // Sort decreasingly diagonal values of s_m & corrensponded eigen vector u_i   


    /* ----- STEP 2 -----*/

    // initialize v 
    make_identity(v, n); 
    
    make_identity(s_n, n); 
    
    // RUN jacobi EVD for A_tA : A_tA = (v)*(s_n)*(v^t) 
    jacobi_evd_n(A_HA, n, v, s_n); 

    eigen_sort(v, s_n, n); // Sort decreasingly diagonal values of s_n & corrensponded eigen vector v_i   

    /* ----- STEP 3 -----*/
    // update Scale Matrix called SIGMA mxn matrix s 
    
    // double s_mxn[m*n];
    make_scale_mxn(s_m, s_n, s, m, n) ;
    std::cout << "------scale matrix SIGMA-----" << std::endl;
    print_mxn_matrix(s, m,n);
    
    /*----- STEP 4 -----*/
    // Jacobi SVD  
    
    // get u_t and v_t 
    // double u_t[m*m] ;
    // double v_t[n*n] ;

    //transpose(u, u_t, m);     
    hermitian_mxn(u, u_H, m, n);
    //transpose(v, v_t, n); 
    hermitian_mxn(v, v_H, m, n);
    
    if (m<=n){
        // (m<n) get U without calculating eig(AAt); 
        // double us[m*n];
        // double newU[m*m];
        
        matmul_mxn(A, v, us, m, n, n, n);
        
        get_u_from_us(us, s, newU, m, n); 
     
        matmul_mxn(newU, s, tmp, m, m, m, n);
        matmul_mxn(tmp, v_H, uSv_H,m,n,n,n); 
   
        for(int i=0; i<m; i++){
            for(int j=0; j<m; j++){
                u[i*m+j] = newU[i*m+j];
            }
        }


    }
    else{ // m> n
        // get V without calculating eig(AtA); 
           
        matmul_mxn(A_H, u, vs_H, n, m, m, m); 
        transpose_mxn(s, s_H, m, n);
              
        get_v_from_vs_t(vs_H, s_H, newV, n, m); 
        
        transpose_mxn(newV, newV_H, n,n); 

        matmul_mxn(u, s, tmp, m, m, m, n);
      
        matmul_mxn(tmp, newV_H, uSv_H,m,n,n,n); 
  
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                v[i*n+j] = newV[i*n+j]; 
            }
        }

        // u is same 
    }

}
