//
// fft.cpp
//
// Code generation for function 'fft'
//

// Include files
#include "fft.h"
#include "FFTImplementationCallback.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>
#include <cstring>

// Function Definitions
namespace coder {
int fft(const creal_T x_data[], int x_size, creal_T y_data[])
{
  array<creal_T, 1U> b_wwc_data;
  array<creal_T, 1U> b_x_data;
  array<creal_T, 1U> b_y;
  array<creal_T, 1U> b_y_data;
  array<creal_T, 1U> r;
  array<creal_T, 1U> y;
  array<double, 2U> costab;
  array<double, 2U> sintab;
  array<double, 2U> sintabinv;
  creal_T wwc_data[2207];
  int nd2;
  int y_size;
  if (x_size == 0) {
    y_size = 0;
  } else {
    double costab1q_data[1105];
    double nt_re;
    int N2blue;
    int costab1q_size_idx_1;
    int i;
    int n2;
    boolean_T useRadix2;
    useRadix2 = ((static_cast<unsigned int>(x_size) &
                  static_cast<unsigned int>(x_size - 1)) == 0U);
    N2blue = internal::fft::FFTImplementationCallback::get_algo_sizes(
        x_size, useRadix2, nd2);
    nt_re = 6.2831853071795862 / static_cast<double>(nd2);
    y_size = static_cast<int>(static_cast<unsigned int>(nd2) >> 1) >> 1;
    costab1q_size_idx_1 = y_size + 1;
    costab1q_data[0] = 1.0;
    nd2 = y_size / 2;
    i = static_cast<unsigned short>(nd2);
    for (int k{0}; k < i; k++) {
      costab1q_data[k + 1] = std::cos(nt_re * (static_cast<double>(k) + 1.0));
    }
    i = nd2 + 1;
    for (int k{i}; k < y_size; k++) {
      costab1q_data[k] = std::sin(nt_re * static_cast<double>(y_size - k));
    }
    costab1q_data[y_size] = 0.0;
    if (!useRadix2) {
      n2 = y_size << 1;
      costab.set_size(1, n2 + 1);
      sintab.set_size(1, n2 + 1);
      costab[0] = 1.0;
      sintab[0] = 0.0;
      sintabinv.set_size(1, n2 + 1);
      for (int k{0}; k < y_size; k++) {
        sintabinv[k + 1] = costab1q_data[(y_size - k) - 1];
      }
      for (int k{costab1q_size_idx_1}; k <= n2; k++) {
        sintabinv[k] = costab1q_data[k - y_size];
      }
      for (int k{0}; k < y_size; k++) {
        costab[k + 1] = costab1q_data[k + 1];
        sintab[k + 1] = -costab1q_data[(y_size - k) - 1];
      }
      for (int k{costab1q_size_idx_1}; k <= n2; k++) {
        costab[k] = -costab1q_data[n2 - k];
        sintab[k] = -costab1q_data[k - y_size];
      }
    } else {
      n2 = y_size << 1;
      costab.set_size(1, n2 + 1);
      sintab.set_size(1, n2 + 1);
      costab[0] = 1.0;
      sintab[0] = 0.0;
      for (int k{0}; k < y_size; k++) {
        costab[k + 1] = costab1q_data[k + 1];
        sintab[k + 1] = -costab1q_data[(y_size - k) - 1];
      }
      for (int k{costab1q_size_idx_1}; k <= n2; k++) {
        costab[k] = -costab1q_data[n2 - k];
        sintab[k] = -costab1q_data[k - y_size];
      }
      sintabinv.set_size(1, 0);
    }
    if (useRadix2) {
      b_x_data.set((creal_T *)&x_data[0], x_size);
      internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
          b_x_data, x_size, costab, sintab, y);
      nd2 = y.size(0);
      y_size = y.size(0);
      for (i = 0; i < nd2; i++) {
        y_data[i] = y[i];
      }
    } else {
      double d;
      double d1;
      double nt_im;
      costab1q_size_idx_1 = (x_size + x_size) - 1;
      n2 = 0;
      wwc_data[x_size - 1].re = 1.0;
      wwc_data[x_size - 1].im = 0.0;
      y_size = x_size << 1;
      for (int k{0}; k <= x_size - 2; k++) {
        nd2 = ((k + 1) << 1) - 1;
        if (y_size - n2 <= nd2) {
          n2 += nd2 - y_size;
        } else {
          n2 += nd2;
        }
        nt_im = -3.1415926535897931 * static_cast<double>(n2) /
                static_cast<double>(x_size);
        i = (x_size - k) - 2;
        wwc_data[i].re = std::cos(nt_im);
        wwc_data[i].im = -std::sin(nt_im);
      }
      i = costab1q_size_idx_1 - 1;
      for (int k{i}; k >= x_size; k--) {
        wwc_data[k] = wwc_data[(costab1q_size_idx_1 - k) - 1];
      }
      y.set_size(x_size);
      y_size = x_size;
      for (i = 0; i < x_size; i++) {
        y_data[i] = y[i];
      }
      for (int k{0}; k < x_size; k++) {
        nd2 = (x_size + k) - 1;
        nt_re = wwc_data[nd2].re;
        nt_im = wwc_data[nd2].im;
        d = x_data[k].im;
        d1 = x_data[k].re;
        y_data[k].re = nt_re * d1 + nt_im * d;
        y_data[k].im = nt_re * d - nt_im * d1;
      }
      i = x_size + 1;
      if (i <= x_size) {
        std::memset(&y_data[i + -1], 0,
                    static_cast<unsigned int>((x_size - i) + 1) *
                        sizeof(creal_T));
      }
      b_y_data.set(&y_data[0], x_size);
      internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
          b_y_data, N2blue, costab, sintab, y);
      b_wwc_data.set(&wwc_data[0], costab1q_size_idx_1);
      internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
          b_wwc_data, N2blue, costab, sintab, r);
      nd2 = y.size(0);
      b_y.set_size(y.size(0));
      for (i = 0; i < nd2; i++) {
        d = y[i].re;
        d1 = r[i].im;
        nt_re = y[i].im;
        nt_im = r[i].re;
        b_y[i].re = d * nt_im - nt_re * d1;
        b_y[i].im = d * d1 + nt_re * nt_im;
      }
      internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
          b_y, N2blue, costab, sintabinv, y);
      if (y.size(0) > 1) {
        nt_re = 1.0 / static_cast<double>(y.size(0));
        nd2 = y.size(0);
        for (i = 0; i < nd2; i++) {
          y[i].re = nt_re * y[i].re;
          y[i].im = nt_re * y[i].im;
        }
      }
      for (int k{x_size}; k <= costab1q_size_idx_1; k++) {
        d = wwc_data[k - 1].re;
        d1 = y[k - 1].im;
        nt_re = wwc_data[k - 1].im;
        nt_im = y[k - 1].re;
        i = k - x_size;
        y_data[i].re = d * nt_im + nt_re * d1;
        y_data[i].im = d * d1 - nt_re * nt_im;
      }
    }
  }
  return y_size;
}

} // namespace coder

// End of code generation (fft.cpp)
