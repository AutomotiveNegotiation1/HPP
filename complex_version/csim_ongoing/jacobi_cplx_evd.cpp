#include <stdio.h>
#include <stdlib.h> 
#include <string.h>
#include <hls_math.h>
#include "jacobi_cplx_evd.h"
#include <iostream>
#include <hls_cordic_apfixed.h>


//double complex_norm(Complex z){
//    return hls::sqrt((z.real() * z.real()) + (z.imag()*z.imag()));
//}

T_48_24 complex_norm(Complex z) {
	T_48_24 real_sq = z.real() * z.real();
	T_48_24 imag_sq = z.imag() * z.imag();
    return hls::sqrt(real_sq + imag_sq);
}

void calculate_sincos(T angle, T_26_2 &sin_val, T_26_2 &cos_val) {
    cordic_apfixed::generic_sincos<48, 24>(angle, sin_val, cos_val);
}

void print_matrix(Complex*A, int N){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
       
            std::cout<< A[i*N+j] << " ";
        }
        std::cout << std::endl; 
    }
    std::cout << std::endl;
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


void jacobi_rotate(Complex *A, Complex *U, int N, int k, int l){
//    double t, c, s;
	T t, c, s;
    Complex c_cplx;
    Complex s_cplx;

    Complex a_kk = A[k*N + k];
    Complex a_ll = A[l*N + l];
    Complex a_kl = A[k*N + l];

//    double phi = arg(a_kl);
    T phi = cordic_apfixed::generic_atan2<48, 24>(a_kl.imag(), a_kl.real());
    T norm = complex_norm(a_kl);
//    T theta = 0.5* cordic_apfixed::generic_atan2(2*complex_norm(a_kl), (a_ll-a_kk).real());
    T theta = T(0.5) * cordic_apfixed::generic_atan2<48, 24>(T(2) * norm, (a_ll - a_kk).real());

    t = cordic_apfixed::generic_tan<48, 24>(theta);
//    c = 1/(hls::sqrt(1+t*t));
    T denom = T(1) + t * t;
    c = T(1) / hls::sqrt(denom);
    s = c*t;

//    c_cplx = Complex(c, 0.0);
    c_cplx = Complex(c, T(0));
//    s_cplx = s* Complex(cos(phi), sin(phi));
//    Complex tmp = Complex(cos(phi), sin(phi));
    ap_fixed<26, 2> sin_phi, cos_phi;
    calculate_sincos(phi, sin_phi, cos_phi);
    Complex tmp = Complex(cos_phi, sin_phi);
    s_cplx = s * tmp;



    for(int h=0; h<N; h++){
        Complex a_hk = A[N*h+k];
        Complex a_hl = A[N*h+l];


        A[N*h+k] = c_cplx*a_hk - x_conj(s_cplx)*a_hl;
        A[N*h+l] = s_cplx*a_hk + c_cplx*a_hl;

        A[N*k+h] = x_conj(A[N*h+k]);
        A[N*l+h] = x_conj(A[N*h+l]);
    }

    A[N*k+l] = A[N*l+k] = Complex(T(0), T(0));

    //A[N*k+k] = c*c*a_kk + s*s*a_ll - 2*s*c*a_kl;
    A[N*k+k] = c_cplx*c_cplx*a_kk + s_cplx*x_conj(s_cplx)*a_ll - T(2)*(s_cplx*c_cplx*x_conj(a_kl)).real();
    A[N*l+l] = s_cplx*x_conj(s_cplx)*a_kk + c_cplx*c_cplx*a_ll + T(2)*(s_cplx*c_cplx*x_conj(a_kl)).real();

    //update U (sample)

    for(int i=0; i<N; i++){

        Complex u_ik = U[i*N + k];
        Complex u_il = U[i*N + l];

        U[i*N+k]= c_cplx*u_ik - x_conj(s_cplx)*u_il;
        U[i*N+l]= s_cplx*u_ik + c_cplx*u_il;
    }


}


// Should Change to get theta in Complex version 
void find_maxDiagOff(Complex* A, int N, int* p, int *q, double* value){

    *value = -1;
    *p=-1; 
    *q=-1; 

    for (int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            
            double tmp = complex_norm(A[i*N+j]); 
            
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
         
//           ret+= complex_norm(A[i*N+j]);
        	ret+= (double)complex_norm(A[i*N+j]);
        }
    }
    return ret; 
}

// Should Change to get theta in Complex version 
void jacobi_evd_m(Complex* A, int N, Complex* u, Complex* s){


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


void jacobi_evd_n(Complex* A, int N, Complex* u, Complex* s){
    
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
