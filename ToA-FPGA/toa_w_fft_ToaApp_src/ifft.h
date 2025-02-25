//
// ifft.h
//
// Code generation for function 'ifft'
//

#ifndef IFFT_H
#define IFFT_H

// Include files
#include "rtwtypes.h"
#include "coder_array.h"
#include <cstddef>
#include <cstdlib>

// Function Declarations
namespace coder {
void ifft(const creal_T x_data[], int x_size, double varargin_1,
          array<creal_T, 1U> &y);

}

#endif
// End of code generation (ifft.h)
