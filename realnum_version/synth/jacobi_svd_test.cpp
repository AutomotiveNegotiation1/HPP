#include "jacobi_svd.h"
#include "svd_test_data.h"
#include <stdbool.h>
#include <hls_stream.h>
#include "ap_axi_sdata.h"

#define m 10
#define n 20

static double usv_t[m*n]; 
static double v_t[n*n];
static double tmp[m*n]; 

typedef hls::axis<double, 0, 0, 1> axis_double;


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
    
    
    transpose_mxn(v, v_t, n,n);
    

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
        return true; 
    }

    else{ // zero_col >0 ; There is no zero _col 
        double u_t[t*t];
        
        double uu_t[t*t];
        
        double u_tu[t*t]; 

        transpose_mxn(u, u_t, t, t);
        matmul_mxn(u, u_t, uu_t, t,t,t,t);
        matmul_mxn(u_t, u, u_tu, t,t,t,t);
        
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

    hls::stream<axis_double> A_in_stream;
    hls::stream<axis_double> U_out_stream;
    hls::stream<axis_double> V_out_stream;
    hls::stream<axis_double> S_out_stream;

    for (int i = 0; i < m * n; i++) {
        axis_double tmp;
        tmp.data = A[i];
        A_in_stream.write(tmp);
    }

    //jacobi_svd(A, u, v, s, m, n);
    jacobi_svd(A_in_stream, U_out_stream, V_out_stream, S_out_stream, m, n);

    static double u[m * m];
    static double v[n * n];
    static double s[m * n];

    for (int i = 0; i < m * m; i++) {
        axis_double tmp = U_out_stream.read();
        u[i] = tmp.data;
    }

    for (int i = 0; i < n * n; i++) {
        axis_double tmp = V_out_stream.read();
        v[i] = tmp.data;
    }

    for (int i = 0; i < m * n; i++) {
        axis_double tmp = S_out_stream.read();
        s[i] = tmp.data;
    }
    
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
    test_data(B); 
    return 0;    
     
}
