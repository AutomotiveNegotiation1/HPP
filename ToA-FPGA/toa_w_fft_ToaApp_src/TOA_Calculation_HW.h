//
// TOA_Calculation_HW.h
//
// Code generation for function 'TOA_Calculation_HW'
//

#ifndef TOA_CALCULATION_HW_H
#define TOA_CALCULATION_HW_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
extern void TOA_Calculation_HW(const creal_T raw_input[17664],
                               const creal_T pskSigOfdm[1024], double simSample,
                               double cp_length, double SC_fil_start,
                               double SC_fil_end,
                               coder::array<creal_T, 3U> &TOA_fft,
                               double *PeakInd_TOA_fft_out, creal_T Azi_data[4],
                               creal_T Ele_data[4]);

#endif
// End of code generation (TOA_Calculation_HW.h)
