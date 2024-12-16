#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "jacobi_evd.h"

#define m 10
#define n 20

void print_matrix(double *A, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%f ", A[i * N + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void make_identity(double *u, int N) {

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                u[i * N + j] = 1;
            } else {
                u[i * N + j] = 0;
            }
        }
    }
}

void transpose(double *A, double *A_t, int N) {

	for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A_t[j * N + i] = A[i * N + j];
        }
    }
}

void matmul(double *A, double *B, double *C, int N) {
#pragma HLS ARRAY_PARTITION variable=A block dim=2
#pragma HLS ARRAY_PARTITION variable=B block dim=2
#pragma HLS ARRAY_PARTITION variable=C block dim=2
    // Assume A, B, C: N by N square matrix
    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1 // 파이프라인 적용 (II=1로 설정하여 각 반복을 파이프라인화)

        for (int j = 0; j < N; j++) {
            C[i * N + j] = 0;
            for (int k = 0; k < N; k++) {
#pragma HLS UNROLL factor=2 // 내부 루프를 4배 언롤하여 병렬화
                C[i * N + j] += A[i * N + k] * B[k * N + j];
            }
        }
    }
}
void jacobi_rotate(double *A, double *U, int N, int k, int l) {
    // 중간 버퍼 선언
    double temp_A_k[20];  // A[k]에 해당하는 행을 임시 저장
    double temp_A_l[20];  // A[l]에 해당하는 행을 임시 저장
    double temp_A_k_update[20];  // A[k] 업데이트용 임시 저장
    double temp_A_l_update[20];  // A[l] 업데이트용 임시 저장
    double temp_U_k[20];  // U[k]에 대한 임시 저장
    double temp_U_l[20];  // U[l]에 대한 임시 저장

    double t, c, s;

    double a_kk = A[k * N + k];
    double a_ll = A[l * N + l];
    double a_kl = A[k * N + l];

    double theta = 0.5 * atan2(2 * a_kl, a_ll - a_kk);
    t = tan(theta);
    c = 1 / (sqrt(1 + t * t));
    s = c * t;
    int h;
    // A의 값들을 중간 버퍼에 먼저 계산
    for (h = 0; h < N; h++) {
        double a_hk = A[N * h + k];
        double a_hl = A[N * h + l];

        // 중간 버퍼에 계산된 값 저장
        temp_A_k[h] = c * a_hk - s * a_hl;
        temp_A_l[h] = c * a_hl + s * a_hk;
    }

    // A[k, l]와 A[l, k]는 0으로 업데이트
    A[N * k + l] = 0;
    A[N * l + k] = 0;

    // A[k, k]와 A[l, l] 업데이트
    A[N * k + k] = c * c * a_kk + s * s * a_ll - 2 * s * c * a_kl;
    A[N * l + l] = s * s * a_kk + c * c * a_ll + 2 * s * c * a_kl;

    // 첫 번째 단계: A[k, h]와 A[l, h]를 업데이트
    // temp_A_k와 temp_A_l을 별도로 업데이트용 중간 버퍼에 저장
    for (h = 0; h < N; h++) {
        temp_A_k_update[h] = temp_A_k[h];  // temp_A_k를 한 번 더 복사
        temp_A_l_update[h] = temp_A_l[h];  // temp_A_l을 한 번 더 복사
    }

    // 첫 번째 루프: A[k, h]와 A[l, h]의 값을 업데이트
    for (h = 0; h < N; h++) {
        A[N * h + k] = temp_A_k_update[h];  // A[k, h] 업데이트
    }
    for (h = 0; h < N; h++) {
		A[N * h + l] = temp_A_l_update[h];  // A[l, h] 업데이트
	}
    // 두 번째 루프: A[k, h]와 A[l, h]의 대칭 행렬을 업데이트
    for (h = 0; h < N; h++) {
        A[N * k + h] = temp_A_k_update[h];  // A[k, h] 대칭 행렬 업데이트
    }

    for (h = 0; h < N; h++) {
            A[N * l + h] = temp_A_l_update[h];  // A[l, h] 대칭 행렬 업데이트
	}
    int i;
    // U의 값들을 중간 버퍼에 먼저 계산
    for (i = 0; i < N; i++) {
        double u_ik = U[i * N + k];
        double u_il = U[i * N + l];

        // 중간 버퍼에 계산된 값 저장
        temp_U_k[i] = c * u_ik - s * u_il;
        temp_U_l[i] = s * u_ik + c * u_il;
    }

    // U 업데이트
    for (i = 0; i < N; i++) {
        U[i * N + k] = temp_U_k[i];
        U[i * N + l] = temp_U_l[i];
    }
}
//void jacobi_rotate(double *A, double *U, int N, int k, int l) {
////#pragma HLS ARRAY_PARTITION variable=A block factor=4 dim=2
////#pragma HLS ARRAY_PARTITION variable=U block factor=4 dim=2
//
//	double t, c, s;
//
//    double a_kk;
//    a_kk = A[k * N + k];
//    double a_ll;
//    a_ll = A[l * N + l];
//    double a_kl;
//	a_kl = A[k * N + l];
//
//    double theta;
//    theta = 0.5 * atan2(2 * a_kl, a_ll - a_kk);
//
//    t = tan(theta);
//    c = 1 / (sqrt(1 + t * t));
//    s = c * t;
//    // Update matrix A
//    for (int h = 0; h < N; h++) {
////#pragma HLS PIPELINE II=1 // II=1로 파이프라인화하여 반복문 간 독립적인 실행 보장
//        double a_hk;
//        a_hk = A[N * h + k];
//        double a_hl;
//        a_hl = A[N * h + l];
//
//        A[N * h + k] = A[N * k + h] = c * a_hk - s * a_hl;
//        A[N * h + l] = A[N * l + h] = c * a_hl + s * a_hk;
//    }
//    A[N * k + l] = A[N * l + k] = 0;
//
//    A[N * k + k] = c * c * a_kk + s * s * a_ll - 2 * s * c * a_kl;
//    A[N * l + l] = s * s * a_kk + c * c * a_ll + 2 * s * c * a_kl;
//
//    // Update U (sample)
//    for (int i = 0; i < N; i++) {
////#pragma HLS PIPELINE II=1 // II=1로 파이프라인화하여 반복문 간 독립적인 실행 보장
//        double u_ik;
//        u_ik = U[i * N + k];
//        double u_il;
//        u_ik = U[i * N + l];
//
//        U[i * N + k] = c * u_ik - s * u_il;
//        U[i * N + l] = s * u_ik + c * u_il;
//    }
//}

void find_maxDiagOff(double *A, int N, int *p, int *q, double *value) {

    *value = -1;
    *p = -1;
    *q = -1;
    for (int i = 0; i < N; i++) {
#pragma HLS PIPELINE II=1 // II=1로 설정하여 파이프라인화
        for (int j = i + 1; j < N; j++) {
            double tmp = A[i * N + j];

            if (tmp < 0) {
                tmp = (-1) * tmp;
            }

            if (*value < tmp) {
                *value = tmp;
                *p = i;
                *q = j;
            }
        }
    }
}

double cal_sumDiagOff(double *A, int N) {

    double ret = 0;
    for (int i = 0; i < N; i++) {

#pragma HLS PIPELINE II=1 // II=1로 설정하여 파이프라인화
        for (int j = i + 1; j < N; j++) {

#pragma HLS UNROLL factor=2 // 내부 루프를 2배 언롤링하여 병렬화
            if (A[i * N + j] < 0.0) {
                ret += (-1) * A[i * N + j];
            } else {
                ret += A[i * N + j];
            }
        }
    }
    return ret;
}

void jacobi_evd(double *A, int N, double *u, double *s) {

    double A_t[N * N];
    double AA_t[N * N];

    transpose(A, A_t, N);
    matmul(A, A_t, AA_t, N);

    int rot_num = 0;
    while (cal_sumDiagOff(A, N) > 1e-10) {
        int p = -1;
        int q = -1;
        double value = -1; // ���� ������ �ʱ�ȭ

        find_maxDiagOff(A, N, &p, &q, &value);
        jacobi_rotate(A, u, N, p, q);

        rot_num++;
    }

    // Update s
    for (int i = 0; i < N; i++) {
        s[i * N + i] = A[i * N + i];
    }
}

void jacobi_evd_m(double *A, int N, double *u, double *s) {
//#pragma HLS INTERFACE ap_stable port=A
//#pragma HLS ARRAY_PARTITION variable=s block dim=1
//#pragma HLS ARRAY_PARTITION variable=A block dim=2
//#pragma HLS ARRAY_PARTITION variable=s block dim=2
//#pragma HLS ARRAY_PARTITION variable=u block dim=1
	static double A_t[m * m];
    static double AA_t[m * m];

    transpose(A, A_t, N);
    matmul(A, A_t, AA_t, N);

    int rot_num = 0;
    static double sumDiagOff;
    sumDiagOff = cal_sumDiagOff(A, N);
    while (sumDiagOff > 1e-10) {
        int p = -1;
        int q = -1;
        double value = -1;

        find_maxDiagOff(A, N, &p, &q, &value);
        jacobi_rotate(A, u, N, p, q);

        rot_num++;
    }

    // Update s
    for (int i = 0; i < N; i++) {

//#pragma HLS UNROLL factor=2  // 대각 성분 업데이트를 병렬화
        s[i * N + i] = A[i * N + i];
    }
}

void jacobi_evd_n(double *A, int N, double *u, double *s) {
//#pragma HLS ARRAY_PARTITION variable=A block dim=2
//#pragma HLS ARRAY_PARTITION variable=s block dim=2
//#pragma HLS ARRAY_PARTITION variable=u block dim=1
	//#pragma HLS INTERFACE ap_stable port=A
//#pragma HLS ARRAY_PARTITION variable=s block dim=1
    static double A_t[n * n];
    static double AA_t[n * n];

    transpose(A, A_t, N);
    matmul(A, A_t, AA_t, N);

    int rot_num = 0;
//#pragma HLS PIPELINE II=1
    static double sumDiagOff;
    sumDiagOff = cal_sumDiagOff(A, N);
    while (sumDiagOff > 1e-10) {
//#pragma HLS UNROLL factor=4
        int p = -1;
        int q = -1;
        double value = -1;

        find_maxDiagOff(A, N, &p, &q, &value);
        jacobi_rotate(A, u, N, p, q);

        rot_num++;
    }

    // Update s
    for (int i = 0; i < N; i++) {
//#pragma HLS UNROLL factor=2  // 대각 성분 업데이트를 병렬화
        s[i * N + i] = A[i * N + i];
    }
}
