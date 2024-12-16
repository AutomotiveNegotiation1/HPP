#include "mat_utils.h"


void print_mxn_matrix(double*A, int m, int n){

    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            printf("%f ",A[i*n+j]); 
        }
        printf("\n");
    }
    printf("\n");
}

void matmul_mxn(double* A, double* B, double*C, int rowA, int colA, int rowB, int colB){

    // A : rowA X colA 
    // B : rowB X colB (colA == rowB)
    // C (= AXB) : rowA X colB 
    if(colA != rowB){
        printf("col# of A != row#B\n"); 
        return; 
    }
    for(int i=0; i<rowA; i++){
        for(int j=0; j<colB; j++){
            double AikBkj = 0 ;
            for(int k=0; k<colA; k++){
                AikBkj += A[i*colA + k] * B[k*colB+ j];
                //printf("i: %d j: %d : %d * %d\n",i,j, A[i*rowA+k], A[k*rowB+j] );
            }
            C[i*colB + j] = AikBkj;
            //printf("C_{%d}_{%d}: %f\n", i, j, AikBkj);
        }
    }
} 

void transpose_mxn(double* A, double* A_t, int m, int n){
    // A: m x n matrix 
    // A_t : n x m matrix
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            A_t[j*m+i] = A[i*n+j]; 
        }
    }
}

void make_scale_mxn(double* s_m, double* s_n, double* s_mxn, int m, int n){


    //int s = -1; 
    if (m>n){
        //s = n;
        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                if(i==j){
                    
                    s_mxn[i*n +j] = sqrt(fabs(s_n[i*n+j])); 
                }
                else{
                    s_mxn[i*n + j] = 0; 
                }
            }
        } 
    }
    else{ // m<=n
        //s = m; 

        for(int i=0; i<m; i++){
            for(int j=0; j<n; j++){
                if(i==j){
                    
                    s_mxn[i*n+j] = sqrt(fabs(s_m[i*m+j])); 
                }
                else{
                    s_mxn[i*n+j]= 0; 
                }
            }
        }
    }

} 

//void swap_vectors(int mth, int nth, double* EV, int N){
//
//    for(int i=0; i<N; i++){
//        // EV[i][col_i] <--> EV[i][col_j]
//        double ev_im = EV[i*N + mth];
//        double ev_in = EV[i*N + nth];
//
//        EV[i*N + mth] = ev_in;
//        EV[i*N + nth] = ev_im;
//    }
//}
void swap_vectors(int mth, int nth, double* EV, int N) {

    // 중간 버퍼 선언
    double temp_mth[10];  // mth 컬럼의 값을 임시 저장할 버퍼
    double temp_nth[10];  // nth 컬럼의 값을 임시 저장할 버퍼

    int i;
    // 첫 번째 단계: mth, nth 인덱스에 해당하는 값을 각각의 임시 버퍼에 저장
    for (i = 0; i < N; i++) {
        temp_mth[i] = EV[i * N + mth];  // mth 컬럼 값 저장
        temp_nth[i] = EV[i * N + nth];  // nth 컬럼 값 저장
    }

    // 두 번째 단계: 임시 버퍼에 저장된 값들을 교환 (각각의 값을 독립적으로 처리)
    for (i = 0; i < N; i++) {
        // 첫 번째로 mth 컬럼에 임시 저장된 nth 값을 교환
        EV[i * N + mth] = temp_nth[i];
    }

    // 이후에 nth 컬럼에 임시 저장된 mth 값을 교환
    for (i = 0; i < N; i++) {
        EV[i * N + nth] = temp_mth[i];
    }
}
void eigen_sort(double* EigenVectors, double* s, int m) {

    // 중간 버퍼 선언
    double s_temp[10 * 10];  // s 배열의 임시 버퍼
    double EigenVectors_temp[10 * 10];  // EigenVectors 배열의 임시 버퍼

    // 원본 배열을 중간 버퍼로 복사
    int i,j;
    for (i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            EigenVectors_temp[i * m + j] = EigenVectors[i * m + j];
            s_temp[i * m + j] = s[i * m + j];
        }
    }

    // 버블 정렬: 중간 버퍼를 사용하여 s와 EigenVectors 배열 수정
    for (i = 0; i < m - 1; i++) {
        for (j = 0; j < m - i - 1; j++) {
            int idx;

            // s[j][j]와 s[j+1][j+1] 값을 중간 버퍼에서 읽어오기
            idx = j * m + j;
            double sjj = s_temp[idx];

            idx = (j + 1) * m + (j + 1);
            double sjj_next = s_temp[idx];

            if (sjj < sjj_next) {

                // 중간 버퍼에서 s 값들을 교환
                idx = (j + 1) * m + (j + 1);
                s_temp[idx] = sjj;
                idx = j * m + j;
                s_temp[idx] = sjj_next;

                // EigenVectors 배열에서의 교환
                swap_vectors(j, j + 1, EigenVectors_temp, m);
            }
        }
    }

    // 최종적으로 중간 버퍼에서 계산된 값을 원본 배열에 반영
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            s[i * m + j] = s_temp[i * m + j];
            EigenVectors[i * m + j] = EigenVectors_temp[i * m + j];
        }
    }
}

void get_u_from_us(double* us, double* s_mxn, double* newU, int m, int n){

    //us : mxn  (us==Av)  
    // s_mxn : mxn 
    // newU : mxm

    if(m<=n){
        for(int j=0; j<m; j++){ // jth column 

#pragma HLS PIPELINE // 파이프라인을 적용하여 성능 향상
            for(int i=0; i<m; i++){ // ith row 
#pragma HLS UNROLL factor=2 // 루프 언롤링을 적용하여 성능 향상
                if (fabs(s_mxn[j*n+j])< 10e-7){
                    s_mxn[j*n+j] = 10e-8; 
                }
                newU[i*m + j] = us[i*n+j] / s_mxn[j*n+j];    
            }
        }
    }    
}


void get_v_from_vs_t(double* vs_t, double* s_mxn_t, double* newV, int n, int m) {
    // 중간 버퍼 정의
    double s_mxn_t_modified[20]; // 각 열에 대한 수정된 값 저장용
//
//    #pragma HLS ARRAY_PARTITION variable=s_mxn_t_modified complete dim=1
//    #pragma HLS ARRAY_PARTITION variable=vs_t complete dim=1
//    #pragma HLS ARRAY_PARTITION variable=newV complete dim=1
    int j;
    if (m > n) {
        // s_mxn_t 배열을 먼저 수정하여 중간 버퍼에 저장
        for (j = 0; j < n; j++) { // j번째 열
            if (fabs(s_mxn_t[j * m + j]) < 10e-7) {
                s_mxn_t_modified[j] = 10e-8; // 중간 버퍼에 수정된 값 저장
            } else {
                s_mxn_t_modified[j] = s_mxn_t[j * m + j]; // 수정되지 않은 값 저장
            }
        }

        // newV 계산
        for (j = 0; j < n; j++) { // j번째 열
//            #pragma HLS PIPELINE
            for (int i = 0; i < n; i++) { // i번째 행
                newV[n * i + j] = vs_t[i * m + j] / s_mxn_t_modified[j]; // 중간 버퍼 사용
            }
        }
    }
}

//void get_v_from_vs_t(double* vs_t, double* s_mxn_t, double* newV, int n, int m){
//    // (v_t)(s_t) : nxm (uA == Sv --> (Sv)^t = (uA)^t)
//    // s_mxn_t : nxm
//    // newV_t : nxn
//
//    if(m>n){
//        for(int j=0; j<n; j++){ //j th column
//
////#pragma HLS PIPELINE // 파이프라인을 적용하여 성능 향상
//            for(int i=0; i<n; i++){ // ith row
//// #pragma HLS UNROLL factor=2 // 루프 언롤링을 적용하여 성능 향상
//                if (fabs(s_mxn_t[j*m+j])< 10e-7){
//                    s_mxn_t[j*m+j] = 10e-8;
//                }
//                newV[n*i + j] = vs_t[i*m+j] / s_mxn_t[j*m+j] ;
//            }
//        }
//    }
//}
