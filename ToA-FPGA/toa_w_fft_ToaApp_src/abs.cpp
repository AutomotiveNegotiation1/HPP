//
// abs.cpp
//
// Code generation for function 'abs'
//

// Include files
#include "abs.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
namespace coder {
void b_abs(const array<creal_T, 1U> &x, array<double, 1U> &y)
{
  int nx_tmp;
  nx_tmp = x.size(0);
  y.set_size(x.size(0));
  for (int k{0}; k < nx_tmp; k++) {
    double a;
    double b;
    a = std::abs(x[k].re);
    b = std::abs(x[k].im);
    if (a < b) {
      a /= b;
      y[k] = b * std::sqrt(a * a + 1.0);
    } else if (a > b) {
      b /= a;
      y[k] = a * std::sqrt(b * b + 1.0);
    } else if (std::isnan(b)) {
      y[k] = rtNaN;
    } else {
      y[k] = a * 1.4142135623730951;
    }
  }
}

} // namespace coder

// End of code generation (abs.cpp)
