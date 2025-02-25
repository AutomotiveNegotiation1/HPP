//
// minOrMax.cpp
//
// Code generation for function 'minOrMax'
//

// Include files
#include "minOrMax.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
namespace coder {
namespace internal {
double maximum(const array<double, 2U> &x)
{
  double ex;
  int last;
  last = x.size(1);
  if (x.size(1) <= 2) {
    if (x.size(1) == 1) {
      ex = x[0];
    } else {
      ex = x[x.size(1) - 1];
      if ((!(x[0] < ex)) && ((!std::isnan(x[0])) || std::isnan(ex))) {
        ex = x[0];
      }
    }
  } else {
    int idx;
    int k;
    if (!std::isnan(x[0])) {
      idx = 1;
    } else {
      boolean_T exitg1;
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last)) {
        if (!std::isnan(x[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (idx == 0) {
      ex = x[0];
    } else {
      ex = x[idx - 1];
      idx++;
      for (k = idx; k <= last; k++) {
        double d;
        d = x[k - 1];
        if (ex < d) {
          ex = d;
        }
      }
    }
  }
  return ex;
}

double maximum(const array<double, 2U> &x, int &idx)
{
  double ex;
  int last_tmp;
  last_tmp = x.size(1);
  if (x.size(1) <= 2) {
    if (x.size(1) == 1) {
      ex = x[0];
      idx = 1;
    } else {
      ex = x[x.size(1) - 1];
      if ((x[0] < ex) || (std::isnan(x[0]) && (!std::isnan(ex)))) {
        idx = x.size(1);
      } else {
        ex = x[0];
        idx = 1;
      }
    }
  } else {
    int k;
    if (!std::isnan(x[0])) {
      idx = 1;
    } else {
      boolean_T exitg1;
      idx = 0;
      k = 2;
      exitg1 = false;
      while ((!exitg1) && (k <= last_tmp)) {
        if (!std::isnan(x[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }
    if (idx == 0) {
      ex = x[0];
      idx = 1;
    } else {
      int i;
      ex = x[idx - 1];
      i = idx + 1;
      for (k = i; k <= last_tmp; k++) {
        double d;
        d = x[k - 1];
        if (ex < d) {
          ex = d;
          idx = k;
        }
      }
    }
  }
  return ex;
}

} // namespace internal
} // namespace coder

// End of code generation (minOrMax.cpp)
