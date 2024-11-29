#include "jacobi_svd.h"
#include "svd_test_data.h"
#include <stdbool.h>
#include "jacobi_config.h"
#include "svd_cplx_test_data.h"

static double usv_t[m*n]; 
static double v_t[n*n];
static double tmp[m*n]; 


bool check_isIdentity(double* u, int t){


    for(int i=0; i<t; i++){
        for(int j=0; j<t; j++){
            
            if(i==j){
                if(fabs(u[t*i+j]-1.0)>1e-5){
                    return false ;
                }
            }
            else{
                if(fabs(u[t*i+j]-0.0)>1e-5){
                    return false; 
                }
            }
            
        }
    }
    return true; 
}

int check_zero_in_diag(double* s){

    int max_num = 0;
    if (m>=n){
        max_num = m;
    }
    else{
        max_num = n;
    }
    //print_mxn_matrix(s, m, n);
    for(int i=0; i<max_num; i++){

        if(fabs(s[i*n+i] - 0.0) < 1e-5){
                return i;
        }
    }
    return -1;
}
bool check_orthogonal_vec(int j1, int j2, double* u, int t){
    // u : mxm  matrix 
    double sum = 0;
    for(int i=0; i<t; i++){
        sum += u[i*t+j1] * u[i*t+j2]; 
    }
    if (fabs(sum - 0.0) > 1e-5){
        return false;  
    }
    else{
        return true ;
    }
}



bool check_svd_result(double* A, double* u, double* v, double* s){
    
    //static double usv_t[m*n]; 
    //static double v_t[n*n];
    
    transpose_mxn(v, v_t, n,n);
    //static double tmp[m*n]; 
    

    matmul_mxn(u, s, tmp, m, m, m, n);
    matmul_mxn(tmp, v_t, usv_t, m,n,n,n); 
    
    // check (u)(s)(v^t) == A 
    for(int i=0; i<m; i++){ 
        for(int j=0 ;j<n; j++){
            if(fabs(usv_t[i*n+j] - A[i*n+j])>1e-7){
                printf("error: %f i: %d j:%d\n", fabs(usv_t[i*n+j] - A[i*n+j]), i, j);
                printf("usv_t[i*n_j]: %f\n", usv_t[i*n+j]);
                printf("A[i*n_j]: %f\n", A[i*n+j]);
                return false;
            }
        }
    }
   
    return true; 

}

bool check_vector_orthogonal(double* u, int t, int val_col_num){
    // check (u)(u^t)= (u^t)(u) = I 
    // check (v)(v^t)= (v^t)(v) = I 

    if(val_col_num != -1){
        for(int j1=0; j1<val_col_num-1; j1++){
            for(int j2=j1+1; j2<val_col_num; j2++){
                if(!check_orthogonal_vec(j1,j2, u, t)){
                    printf("orthogonal vector : False\n"); 
                    return false;
                }; 
            }
        }
        //printf("chk orthogonal mat: true\n"); 
        return true; 
    }

    else{ // zero_col >0 ; There is no zero _col 
        double u_t[t*t];
        
        double uu_t[t*t];
        
        double u_tu[t*t]; 

        transpose_mxn(u, u_t, t, t);
        matmul_mxn(u, u_t, uu_t, t,t,t,t);
        matmul_mxn(u_t, u, u_tu, t,t,t,t);
        
        /*
        printf("----check u or v \n");
        print_matrix(u, m);

        printf("----check Identity uut or vvt\n");
        print_matrix(uu_t, m);
        printf("----check Identity utu or vvt\n");
        print_matrix(u_tu, m);
        */
        if(check_isIdentity(uu_t, t) && check_isIdentity(u_tu, t)){
            return true; 
        }
        else{
            return false; 
        }
        
    }
    
}


bool check_vector_orthogonal_m(double* u, int val_col_num){
    // check (u)(u^t)= (u^t)(u) = I 
    // check (v)(v^t)= (v^t)(v) = I 

    if(val_col_num != -1){
        for(int j1=0; j1<val_col_num-1; j1++){
            for(int j2=j1+1; j2<val_col_num; j2++){
                if(!check_orthogonal_vec(j1,j2, u, m)){
                    printf("orthogonal vector : False\n"); 
                    return false;
                }; 
            }
        }
        //printf("chk orthogonal mat: true\n"); 
        return true; 
    }

    else{ // zero_col >0 ; There is no zero _col 
        static double u_t[m*m];
        
        static double uu_t[m*m];
        
        static double u_tu[m*m]; 

        transpose_mxn(u, u_t, m, m);
        matmul_mxn(u, u_t, uu_t, m,m,m,m);
        matmul_mxn(u_t, u, u_tu, m,m,m,m);
        
        if(check_isIdentity(uu_t, m) && check_isIdentity(u_tu, m)){
            return true; 
        }
        else{
            return false; 
        }
        
    }
    
}


bool check_vector_orthogonal_n(double* v, int val_col_num){
    // check (u)(u^t)= (u^t)(u) = I 
    // check (v)(v^t)= (v^t)(v) = I 

    if(val_col_num != -1){
        for(int j1=0; j1<val_col_num-1; j1++){
            for(int j2=j1+1; j2<val_col_num; j2++){
                if(!check_orthogonal_vec(j1,j2, v, n)){
                    printf("orthogonal vector : False\n"); 
                    return false;
                }; 
            }
        }
        //printf("chk orthogonal mat: true\n"); 
        return true; 
    }

    else{ // zero_col >0 ; There is no zero _col 
        static double v_t[n*n];
        
        static double vv_t[n*n];
        
        static double v_tv[n*n]; 

        transpose_mxn(v, v_t, n, n);
        matmul_mxn(v, v_t, vv_t, n,n,n,n);
        matmul_mxn(v_t, v, v_tv, n,n,n,n);
        

        if(check_isIdentity(vv_t, n) && check_isIdentity(v_tv, n)){
            return true; 
        }
        else{
            return false; 
        }
        
    }
    
}

void test_data(double* A){

    static double u[m*m];
    static double v[n*n];
    static double s[m*n];

    //jacobi_svd(A, u, v, s, m, n);
    jacobi_svd(A, u, v, s, m, n);
    
    printf("-----u ------\n");
    print_matrix(u, m);
    
    printf("-----v------\n");
    print_matrix(v, n);
    

    
    if(check_svd_result(A, u, v, s)){
        int valid_col_num = check_zero_in_diag(s);

        if(check_vector_orthogonal_m(u, valid_col_num) && check_vector_orthogonal_n(v, valid_col_num)){
            printf("Test  was successfully completed.\n"); 
        };

    }
    else{
        printf("Test  had an error.\n");
    }
    
}



int main(){
    //test_data(A); 
    test_data(B); 
    //test_data(C, 20, 10); 
    //test_data(D, 5, 5); 
    //test_data(E, 30, 30); 
    /*
    test_data(F, 30, 50); 
    test_data(G, 50, 30); 
    test_data(H, 40, 40); 
    test_data(J, 60, 40); 
    test_data(K, 70, 70); 
    test_data(L, 5, 5); 
    */


    /*
    // #1: Test Data : A 
    int m1 = 5; 
    int n1 = 5; 
    double u1[m1*m1];
    double v1[n1*n1];
    double s1[m1*n1];
    jacobi_svd(A, u1, v1, s1, m1, n1);

    //printf("-----u ------\n");
    //print_matrix(u1, m1);
    if(check_svd_result(A, u1, v1, s1, m1, n1)){
        printf("Test 1 was successfully completed.\n"); 
    }
    else{
        printf("Test 1 had an error.\n");
    }
    */

    return 0;    
     
}
