#include "jacobi_cplx_svd.h"
#include <stdbool.h>
#include "jacobi_config.h"
#include "svd_cplx_test_data.h"



bool check_svd_result(Complex* A, Complex* u, Complex* v, Complex* s){
    
    static Complex usv_H[m*n]; 
    static Complex v_H[n*n];
    static Complex tmp[m*n]; 
    
    hermitian_mxn(v, v_H, n,n);
    

    matmul_mxn(u, s, tmp, m, m, m, n);
    matmul_mxn(tmp, v_H, usv_H, m,n,n,n); 
    
    std::cout << "-----Jacobi SVD result: usv_H-----" << std::endl;
    print_mxn_matrix(usv_H, m,n);
    
    // check (u)(s)(v^t) == A 
    for(int i=0; i<m; i++){ 
        for(int j=0 ;j<n; j++){
            if(complex_norm(usv_H[i*n+j] - A[i*n+j])>1e-7){
                //printf("error: %f i: %d j:%d\n", fabs(usv_H[i*n+j] - A[i*n+j]), i, j);
                //printf("usv_t[i*n_j]: %f\n", usv_H[i*n+j]);
                //printf("A[i*n_j]: %f\n", A[i*n+j]);
                return false;
            }
        }
    }
   
    return true; 

}
void test_data(Complex* A){

    static Complex u[m*m];
    static Complex v[n*n];
    static Complex s[m*n];

    //jacobi_svd(A, u, v, s, m, n);
    jacobi_svd(A, u, v, s);
    
    std::cout << "-----u ------" << std::endl;
    print_matrix(u, m);
    
    std::cout << "-----v------" << std::endl;
    print_matrix(v, n);
    
    if(check_svd_result(A, u, v, s)){

        std::cout << "Test  was successfully completed." << std::endl;
    }
    
    else{
        std::cout << "Test  had an error." << std::endl; 
    }
    
}



int main(){
   
    //test_data(cA); 
    test_data(cB);
   
    return 0;    
     
}
