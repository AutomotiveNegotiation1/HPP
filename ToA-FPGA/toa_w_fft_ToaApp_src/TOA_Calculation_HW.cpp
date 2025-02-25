//
// TOA_Calculation_HW.cpp
//
// Code generation for function 'TOA_Calculation_HW'
//

// Include files
#include "TOA_Calculation_HW.h"
#include "TOA_Calculation_HW_data.h"
#include "abs.h"
#include "any1.h"
#include "fft.h"
#include "findpeaks.h"
#include "ifft.h"
#include "minOrMax.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <stdint.h>
#include "xaxidma.h"
#include "xparameters.h"
#include "xil_printf.h"
#include "xaxidma_hw.h"
#include "xil_io.h"
#include "sleep.h"
#include "xtime_l.h"
#include <iostream>


int init_dma(XAxiDma *dma_instance, u16 device_id) {
    XAxiDma_Config *dma_config;
    dma_config = XAxiDma_LookupConfig(device_id);
    if (!dma_config) {
        xil_printf("No DMA config found for Device ID %d.\n", device_id);
        return XST_FAILURE;
    }
    if (XAxiDma_CfgInitialize(dma_instance, dma_config) != XST_SUCCESS) {
        xil_printf("DMA initialization failed for Device ID %d.\n", device_id);
        return XST_FAILURE;
    }
    if (XAxiDma_HasSg(dma_instance)) {
        xil_printf("Device ID %d configured as SG mode. SG mode is not supported.\n", device_id);
        return XST_FAILURE;
    }
    xil_printf("DMA initialized successfully for Device ID %d.\n", device_id);
    return XST_SUCCESS;
}


// Function Definitions
void TOA_Calculation_HW(const creal_T raw_input[17664],
                        const creal_T pskSigOfdm[1024], double simSample,
                        double cp_length, double SC_fil_start,
                        double SC_fil_end, coder::array<creal_T, 3U> &TOA_fft,
                        double *PeakInd_TOA_fft_out, creal_T Azi_data[4],
                        creal_T Ele_data[4])
{
  coder::array<creal_T, 1U> b_TOA_fft;
  coder::array<double, 2U> TOA_profile;
  coder::array<double, 2U> b_TOA_profile;
  coder::array<double, 2U> peak_inds;
  coder::array<double, 2U> peak_vals;
  coder::array<double, 1U> r;
  coder::array<int, 2U> r2;
  coder::array<boolean_T, 2U> valid_mask;
  creal_T rxSigFreq_data[1104];
  double ai;
  double peak_index;
  double y;
  int b_loop_ub_tmp;
  int i;
  int i1;
  int i2;
  int i3;
  int loop_ub_tmp;
  int pskSigOfdm_re_tmp;
  int rxSigFreq_size;
  int unnamed_idx_1;


  XTime start_time, end_time1, end_time2, start_time_i, end_time_i;
  XAxiDma AxiDma_FFT_in,AxiDma_FFT_out;
  int status, status1, status2;

  xil_printf("Initializing DMA for FFT...\n");
  status = init_dma(&AxiDma_FFT_in, XPAR_AXI_DMA_3_DEVICE_ID);
  if (status != XST_SUCCESS) {
      xil_printf("Failed to initialize DMA for FFT in.\n");
      return;
  }
  status = init_dma(&AxiDma_FFT_out, XPAR_AXI_DMA_2_DEVICE_ID);
  if (status != XST_SUCCESS) {
      xil_printf("Failed to initialize DMA for FFT out.\n");
      return;
  }

  XAxiDma_IntrDisable(&AxiDma_FFT_in, XAXIDMA_IRQ_ALL_MASK, XAXIDMA_DMA_TO_DEVICE);
  XAxiDma_IntrDisable(&AxiDma_FFT_out, XAXIDMA_IRQ_ALL_MASK, XAXIDMA_DEVICE_TO_DMA);


  //  입력:
  //    raw_input    : 수신된 신호 (CP 포함), 크기: (simSample+cp_length) x
  //    RX_Azi x RX_Ele pskSigOfdm   : 전송된 OFDM 심볼 (주파수 도메인), 길이
  //    simSample인 벡터 또는 스칼라 simSample    : OFDM 심볼의 FFT/IFFT 포인트
  //    수 (CP 제외) cp_length    : 사이클릭 프리픽스 길이 SC_fil_start : 유효
  //    subcarrier 시작 인덱스 (정수, ≥1) SC_fil_end   : 유효 subcarrier 끝
  //    인덱스 (정수, ≤ simSample)
  //
  //  출력:
  //    TOA_fft             : 시간 도메인 TOA 프로파일, 크기: simSample x RX_Azi
  //    x RX_Ele PeakInd_TOA_fft_out : 보정된 피크 인덱스 (오프셋 및 circular
  //    보정 적용) Azi_data            : 추출된 azimuth DOA 데이터 (벡터)
  //    Ele_data            : 추출된 elevation DOA 데이터 (벡터)
  //  사전 처리
  //  TOA 프로파일 계산 (각 안테나별 처리)
  loop_ub_tmp = static_cast<int>(simSample);
  TOA_fft.set_size(loop_ub_tmp, 4, 4);
  b_loop_ub_tmp = static_cast<int>(simSample) << 4;
  for (i = 0; i < b_loop_ub_tmp; i++) {
    TOA_fft[i].re = 0.0;
    TOA_fft[i].im = 0.0;
  }
  if (cp_length + 1.0 > 1104.0) {
    i = 0;
    i1 = 0;
  } else {
    i = static_cast<int>(cp_length + 1.0) - 1;
    i1 = 1104;
  }
  b_loop_ub_tmp = i1 - i;
  y = std::sqrt(simSample);
  if (SC_fil_start > SC_fil_end) {
    i1 = 0;
    unnamed_idx_1 = 0;
    i2 = 0;
    i3 = 0;
  } else {
    i1 = static_cast<int>(SC_fil_start) - 1;
    unnamed_idx_1 = static_cast<int>(SC_fil_end);
    i2 = static_cast<int>(SC_fil_start) - 1;
    i3 = static_cast<int>(SC_fil_start) - 1;
  }
  unnamed_idx_1 -= i1;

  XTime_GetTime(&start_time_i);
  for (int arr_Azi_i{0}; arr_Azi_i < 4; arr_Azi_i++) {
    for (int arr_ELe_i{0}; arr_ELe_i < 4; arr_ELe_i++) {
      creal_T raw_input_data[1104];
      double re;
      //  1. CP 제거
      //  2. FFT 수행 (정규화 포함)
      for (i1 = 0; i1 < b_loop_ub_tmp; i1++) {
        raw_input_data[i1] =
            raw_input[((i + i1) + 1104 * arr_Azi_i) + 4416 * arr_ELe_i];
      }

      std::complex<float> dma_input[1104];
      std::complex<float> dma_output[1024];

      XTime_GetTime(&start_time);

      for (int i = 0; i < 1104; i++) {
          dma_input[i] = std::complex<float>(raw_input_data[i].re, raw_input_data[i].im);
      }

      xil_printf("Sending input data to FFT IP via DMA...\n");

      Xil_DCacheFlushRange((UINTPTR)dma_input, sizeof(std::complex<float>) * 1024);
      Xil_DCacheFlushRange((UINTPTR)dma_output, sizeof(std::complex<float>) * 1024);

      status1 = XAxiDma_SimpleTransfer(&AxiDma_FFT_out, (UINTPTR)dma_output, sizeof(std::complex<float>) * 1024, XAXIDMA_DEVICE_TO_DMA);
          if (status1 != XST_SUCCESS) {
              xil_printf("DMA transfer for FFT output failed.\n");
              return;
          }

      status2 = XAxiDma_SimpleTransfer(&AxiDma_FFT_in, (UINTPTR)dma_input, sizeof(std::complex<float>) * 1024, XAXIDMA_DMA_TO_DEVICE);
      if (status2 != XST_SUCCESS) {
          xil_printf("DMA transfer for FFT input failed.\n");
          return;
      }

      while (XAxiDma_Busy(&AxiDma_FFT_in, XAXIDMA_DMA_TO_DEVICE)) {
          xil_printf("Waiting for FFT input DMA transfer to complete...\n");
      }

      while (XAxiDma_Busy(&AxiDma_FFT_out, XAXIDMA_DEVICE_TO_DMA)) {
          xil_printf("Waiting for FFT output DMA transfer to complete...\n");
      }

      Xil_DCacheInvalidateRange((UINTPTR)dma_output, sizeof(std::complex<float>) * 1024);

      for (int i = 0; i < 1024; i++) {
    	  rxSigFreq_data[i].re = dma_output[i].real();
    	  rxSigFreq_data[i].im = dma_output[i].imag();
      }

      XTime_GetTime(&end_time1);

      double elapsed_time = 1.0 * (end_time1 - start_time) / COUNTS_PER_SECOND;
      printf("FFT Elapsed time: %f seconds\n", elapsed_time);


//      rxSigFreq_size =
//          coder::fft(raw_input_data, b_loop_ub_tmp, rxSigFreq_data);
      rxSigFreq_size = 1024;

      for (i1 = 0; i1 < rxSigFreq_size; i1++) {
        peak_index = rxSigFreq_data[i1].re;
        ai = rxSigFreq_data[i1].im;
        if (ai == 0.0) {
          re = peak_index / y;
          peak_index = 0.0;
        } else if (peak_index == 0.0) {
          re = 0.0;
          peak_index = ai / y;
        } else {
          re = peak_index / y;
          peak_index = ai / y;
        }
        rxSigFreq_data[i1].re = re;
        rxSigFreq_data[i1].im = peak_index;
      }
      //  3. 매치드 필터링: 관심 subcarrier 영역에 대해 전송 신호의 켤레 곱함
      for (i1 = 0; i1 < unnamed_idx_1; i1++) {
        double d;
        pskSigOfdm_re_tmp = i2 + i1;
        peak_index = pskSigOfdm[pskSigOfdm_re_tmp].re;
        re = -pskSigOfdm[pskSigOfdm_re_tmp].im;
        ai = rxSigFreq_data[pskSigOfdm_re_tmp].re;
        d = rxSigFreq_data[pskSigOfdm_re_tmp].im;
        raw_input_data[i1].re = ai * peak_index - d * re;
        raw_input_data[i1].im = ai * re + d * peak_index;
      }
      for (i1 = 0; i1 < unnamed_idx_1; i1++) {
        rxSigFreq_data[i3 + i1] = raw_input_data[i1];
      }
      //  4. IFFT 수행하여 시간 도메인 TOA 프로파일 생성
      coder::ifft(rxSigFreq_data, rxSigFreq_size, simSample, b_TOA_fft);
      for (i1 = 0; i1 < loop_ub_tmp; i1++) {
        TOA_fft[(i1 + TOA_fft.size(0) * arr_Azi_i) +
                TOA_fft.size(0) * 4 * arr_ELe_i] = b_TOA_fft[i1];
      }
    }
  }
  XTime_GetTime(&end_time_i);

  double elapsed_time_i = 1.0 * (end_time_i - start_time_i) / COUNTS_PER_SECOND;
  printf("block Elapsed time: %f seconds\n", elapsed_time_i);
  //  피크 검출: 오프셋 및 circular 보정 적용
  //  첫 번째 안테나 (1,1)의 TOA 프로파일 사용
  b_TOA_fft.set_size(loop_ub_tmp);
  for (i = 0; i < loop_ub_tmp; i++) {
    b_TOA_fft[i] = TOA_fft[i];
  }
  coder::b_abs(b_TOA_fft, r);
  pskSigOfdm_re_tmp = r.size(0);
  TOA_profile.set_size(1, r.size(0));
  for (i = 0; i < pskSigOfdm_re_tmp; i++) {
    TOA_profile[i] = r[i];
  }
  //  1 x simSample 벡터
  peak_index = coder::internal::maximum(TOA_profile);
  TOA_profile.set_size(1, TOA_profile.size(1));
  pskSigOfdm_re_tmp = TOA_profile.size(1) - 1;
  unnamed_idx_1 = (TOA_profile.size(1) / 2) << 1;
  rxSigFreq_size = unnamed_idx_1 - 2;
  for (i = 0; i <= rxSigFreq_size; i++) {
      TOA_profile[i] = TOA_profile[i] / peak_index;
  }
//  for (i = 0; i <= rxSigFreq_size; i += 2) {
//    __m128d r1;
//    r1 = _mm_loadu_pd(&TOA_profile[i]);
//    _mm_storeu_pd(&TOA_profile[i], _mm_div_pd(r1, _mm_set1_pd(peak_index)));
//  }
  for (i = unnamed_idx_1; i <= pskSigOfdm_re_tmp; i++) {
    TOA_profile[i] = TOA_profile[i] / peak_index;
  }
  //  정규화
  //  --- 오프셋 보정을 위한 확장 ---
  //  원본 코드에서는 TOA 프로파일의 끝 50개 샘플을 앞에 붙여 확장한 후 피크
  //  검출을 수행함. 보정에 사용할 샘플 수 (필요에 따라 조정)
  if (TOA_profile.size(1) - 49 > TOA_profile.size(1)) {
    i = -1;
    i1 = -1;
  } else {
    i = TOA_profile.size(1) - 51;
    i1 = TOA_profile.size(1) - 1;
  }
  //  extended_profile의 길이는 extension_length + simSample
  //  피크 검출 (findpeaks)
  pskSigOfdm_re_tmp = i1 - i;
  b_TOA_profile.set_size(1, pskSigOfdm_re_tmp + TOA_profile.size(1));
  for (unnamed_idx_1 = 0; unnamed_idx_1 < pskSigOfdm_re_tmp; unnamed_idx_1++) {
    b_TOA_profile[unnamed_idx_1] = TOA_profile[(i + unnamed_idx_1) + 1];
  }
  pskSigOfdm_re_tmp = TOA_profile.size(1);
  for (unnamed_idx_1 = 0; unnamed_idx_1 < pskSigOfdm_re_tmp; unnamed_idx_1++) {
    b_TOA_profile[(unnamed_idx_1 + i1) - i] = TOA_profile[unnamed_idx_1];
  }
  coder::findpeaks(b_TOA_profile, peak_vals, peak_inds);
  //  만약 피크가 검출되지 않으면 fallback: 원래 프로파일에서 최대값 선택
  if (peak_inds.size(1) == 0) {
    coder::internal::maximum(TOA_profile, unnamed_idx_1);
    peak_index = unnamed_idx_1;
  } else {
    //  --- 후보 피크들의 실제 인덱스 계산 ---
    //  extended_profile의 인덱스 1:extension_length는 원래 TOA_profile의
    //  인덱스 (simSample - extension_length + 1) ~ simSample에 해당하고,
    //  인덱스 (extension_length+1) ~ (extension_length+simSample)는 원래 인덱스
    //  1 ~ simSample에 해당함. -1을 왜 추가로 해주는진 모르겠음...;;
    TOA_profile.set_size(1, peak_inds.size(1));
    loop_ub_tmp = peak_inds.size(1);
    for (i = 0; i < loop_ub_tmp; i++) {
      TOA_profile[i] = 0.0;
    }
    i = peak_inds.size(1);
    for (unnamed_idx_1 = 0; unnamed_idx_1 < i; unnamed_idx_1++) {
      ai = peak_inds[unnamed_idx_1];
      if (ai > 50.0) {
        TOA_profile[unnamed_idx_1] = (ai - 50.0) - 1.0;
      } else {
        TOA_profile[unnamed_idx_1] = ((ai + simSample) - 50.0) - 1.0;
      }
    }
    //  delay는 0 이상 simSample/2 이하 범위 내여야 한다고 가정
    valid_mask.set_size(1, peak_inds.size(1));
    peak_index = simSample / 2.0;
    for (i = 0; i < loop_ub_tmp; i++) {
      valid_mask[i] = (TOA_profile[i] <= peak_index);
    }
    if (coder::any(valid_mask)) {
      //  유효 후보들 중 진폭이 가장 큰 것을 선택
      unnamed_idx_1 = 0;
      for (rxSigFreq_size = 0; rxSigFreq_size < loop_ub_tmp; rxSigFreq_size++) {
        if (valid_mask[rxSigFreq_size]) {
          unnamed_idx_1++;
        }
      }
      r2.set_size(1, unnamed_idx_1);
      unnamed_idx_1 = 0;
      for (rxSigFreq_size = 0; rxSigFreq_size < loop_ub_tmp; rxSigFreq_size++) {
        if (valid_mask[rxSigFreq_size]) {
          r2[unnamed_idx_1] = rxSigFreq_size;
          unnamed_idx_1++;
        }
      }
      pskSigOfdm_re_tmp = r2.size(1);
      b_TOA_profile.set_size(1, r2.size(1));
      for (i = 0; i < pskSigOfdm_re_tmp; i++) {
        b_TOA_profile[i] = peak_vals[r2[i]];
      }
      coder::internal::maximum(b_TOA_profile, unnamed_idx_1);
      peak_index = TOA_profile[r2[unnamed_idx_1 - 1]];
    } else {
      //  유효 후보가 없으면 전체 중 최대값 선택 후 circular 보정
      coder::internal::maximum(peak_vals, unnamed_idx_1);
      peak_index = TOA_profile[unnamed_idx_1 - 1];
      if (peak_index > simSample / 2.0) {
        peak_index -= simSample;
      }
    }
  }
  //  DOA 데이터 추출 (피크 인덱스에서 azimuth 및 elevation 데이터 선택)
  Azi_data[0] = TOA_fft[static_cast<int>(peak_index) - 1];
  Ele_data[0] = TOA_fft[static_cast<int>(peak_index) - 1];
  Azi_data[1] = TOA_fft[(static_cast<int>(peak_index) + TOA_fft.size(0)) - 1];
  Ele_data[1] =
      TOA_fft[(static_cast<int>(peak_index) + TOA_fft.size(0) * 4) - 1];
  Azi_data[2] =
      TOA_fft[(static_cast<int>(peak_index) + TOA_fft.size(0) * 2) - 1];
  Ele_data[2] =
      TOA_fft[(static_cast<int>(peak_index) + TOA_fft.size(0) * 4 * 2) - 1];
  Azi_data[3] =
      TOA_fft[(static_cast<int>(peak_index) + TOA_fft.size(0) * 3) - 1];
  Ele_data[3] =
      TOA_fft[(static_cast<int>(peak_index) + TOA_fft.size(0) * 4 * 3) - 1];
  *PeakInd_TOA_fft_out = peak_index;

  XTime_GetTime(&end_time2);

  double elapsed_time = 1.0 * (end_time2 - start_time_i) / COUNTS_PER_SECOND;
  printf("TOA Elapsed time: %f seconds\n", elapsed_time);
}

// End of code generation (TOA_Calculation_HW.cpp)
