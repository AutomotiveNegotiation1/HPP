fft_ide : fft만 ide에서 사용하는 예제용 SDk IDE 코드 
fft_vivado : FFT가 포함된 SDk IDE용 Vivado 코드 
toa_w_fft : fft 블락을 포함하는 SDK IDE 코드 
  일반 C ToA 변경법 : toa_w_fft/ToaApp/src/ToA_Calculation_HW.cpp 내의 
    - line155:205 주석처리 (FFT 부분, for (i1=0;i<b_loop_up_bm;i1++) { ~~printf("FFT Elapsed time: "); 
    - line208:209 주석해제 (rxSigFreq_size = coder::fft(raw_input_data, b_loop_ub_tmp, rxSigFreq_data); 
                        
FFT 블럭 사용 시간 (ZCU102, 1.2GHz) 
  * 4x4 루프를 병렬화하지 않았음 
  FFT Elapsed time: 0.005347 seconds 
  block Elapsed time: 0.175689 seconds 
  ToA Elapsed time: 0.1777384 seconds 
  Total Elapsed time: 0.216471 seconds 

FFT 블럭 사용 안하고 PL 사용 (ZCU102, 1.2GHz) 
  block Elapsed time: 0.074915 seconds 
  ToA Elapsed time: 0.075830 seconds 
  Total Elapsed time: 0.112563 seconds
*전체 ToA 시간 중 FFT가 99%를 차지함. 이에 HW가속화 필요. 
* 1104x1x1에 대한 안테나 array 개별에 대해서 fft 수행하여 4x4 2중 루프로 16번 수행함
