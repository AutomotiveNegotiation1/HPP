// Include files
#include "main.h"
#include "TOA_Calculation_HW.h"
#include "TOA_Calculation_HW_terminate.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <iostream>
#include "mat_data.h"
#include <fstream>
#include "xtime_l.h"

// Function Declarations
static void argInit_1024x1_creal_T(creal_T result[1024]);

static void argInit_1104x4x4_creal_T(creal_T result[17664]);

static creal_T argInit_creal_T();

static double argInit_real_T();

// Function Definitions
static void argInit_1024x1_creal_T(creal_T result[1024])
{
  // Loop over the array to initialize each element.
  for (int idx0{0}; idx0 < 1024; idx0++) {
    // Set the value of the array element.
    // Change this value to the value that the application requires.
    result[idx0] = argInit_creal_T();
  }
}

static void argInit_1104x4x4_creal_T(creal_T result[17664])
{
  // Loop over the array to initialize each element.
  for (int idx0{0}; idx0 < 1104; idx0++) {
    for (int idx1{0}; idx1 < 4; idx1++) {
      for (int idx2{0}; idx2 < 4; idx2++) {
        // Set the value of the array element.
        // Change this value to the value that the application requires.
        result[(idx0 + 1104 * idx1) + 4416 * idx2] = argInit_creal_T();
      }
    }
  }
}

static creal_T argInit_creal_T()
{
  creal_T result;
  double re_tmp;
  // Set the value of the complex variable.
  // Change this value to the value that the application requires.
  re_tmp = argInit_real_T();
  result.re = re_tmp;
  result.im = re_tmp;
  return result;
}

static double argInit_real_T()
{
  return 0.0;
}


int main(int, char **)
{
	XTime start_time_tot, end_time_tot;

  // The initialize function is being called automatically from your entry-point
  // function. So, a call to initialize is not included here. Invoke the
  // entry-point functions.
  // You can call entry-point functions multiple times.
  std::cout << "Main function started." << std::endl;
  XTime_GetTime(&start_time_tot);
  main_TOA_Calculation_HW();
  // Terminate the application.
  // You do not need to do this more than one time.
  TOA_Calculation_HW_terminate();
  XTime_GetTime(&end_time_tot);

  double elapsed_time = 1.0 * (end_time_tot - start_time_tot) / COUNTS_PER_SECOND;
  printf("Total Elapsed time: %f seconds\n", elapsed_time);

  return 0;
}


// void main_TOA_Calculation_HW()
// {
//     static creal_T dcv[17664];  // raw_input
//     coder::array<creal_T, 3U> TOA_fft;
//     creal_T dcv1[1024];  // pskSigOfdm
//     creal_T Azi_data[4];
//     creal_T Ele_data[4];

//     double simSample, cp_length, SC_fil_start, SC_fil_end;
//     double PeakInd_TOA_fft_out;


//     std::cout << "main_TOA_Calculation_HW Call." << std::endl;
//     // load data from `.mat` file
//     loadMatFile("sim_input_data.mat", dcv, dcv1, simSample, cp_length, SC_fil_start, SC_fil_end);

//     TOA_Calculation_HW(dcv, dcv1, simSample, cp_length, SC_fil_start, SC_fil_end, 
//                        TOA_fft, &PeakInd_TOA_fft_out, Azi_data, Ele_data);

//     std::cout << "TOA_Calculation_HW execution completed." << std::endl;

//     std::cout << "PeakInd_TOA_fft_out: " << PeakInd_TOA_fft_out << std::endl;
//     std::cout << "Azi_data: ";
//     for (int i = 0; i < 4; i++) {
//         std::cout << "(" << Azi_data[i].re << ", " << Azi_data[i].im << ") ";
//     }
//     std::cout << std::endl;
//     std::cout << "Ele_data: ";
//     for (int i = 0; i < 4; i++) {
//         std::cout << "(" << Ele_data[i].re << ", " << Ele_data[i].im << ") ";
//     }
//     std::cout << std::endl;
// }


void main_TOA_Calculation_HW()
{
    coder::array<creal_T, 3> TOA_fft;  
    creal_T Azi_data[4];
    creal_T Ele_data[4];
    double PeakInd_TOA_fft_out;

    std::cout << "main_TOA_Calculation_HW Call." << std::endl;

    TOA_Calculation_HW(dcv, dcv1, simSample, cp_length, SC_fil_start, SC_fil_end, 
                       TOA_fft, &PeakInd_TOA_fft_out, Azi_data, Ele_data);

     std::cout << "TOA_fft size: [" << TOA_fft.size(0) << ", " 
               << TOA_fft.size(1) << ", " << TOA_fft.size(2) << "]" << std::endl;

    std::cout << "TOA_Calculation_HW execution completed." << std::endl;
    std::cout << "PeakInd_TOA_fft_out: " << PeakInd_TOA_fft_out << std::endl;
    
    std::cout << "Azi_data: ";
    for (int i = 0; i < 4; i++) {
        std::cout << "(" << Azi_data[i].re << ", " << Azi_data[i].im << ") ";
    }
    std::cout << std::endl;
    
    std::cout << "Ele_data: ";
    for (int i = 0; i < 4; i++) {
        std::cout << "(" << Ele_data[i].re << ", " << Ele_data[i].im << ") ";
    }
    std::cout << std::endl;
}

