//
// ifft.cpp
//
// Code generation for function 'ifft'
//

// Include files
#include "ifft.h"
#include "FFTImplementationCallback.h"
#include "TOA_Calculation_HW_data.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
namespace coder {
void ifft(const creal_T x_data[], int x_size, double varargin_1,
          array<creal_T, 1U> &y)
{
  array<creal_T, 1U> b_fv;
  array<creal_T, 1U> b_x_data;
  array<creal_T, 1U> fv;
  array<creal_T, 1U> wwc;
  array<double, 2U> costab;
  array<double, 2U> costab1q;
  array<double, 2U> sintab;
  array<double, 2U> sintabinv;
  int nd2;
  int nfft_tmp_tmp;
  nfft_tmp_tmp = static_cast<int>(varargin_1);
  if ((x_size == 0) || (static_cast<int>(varargin_1) == 0)) {
    y.set_size(nfft_tmp_tmp);
    for (int i{0}; i < nfft_tmp_tmp; i++) {
      y[i].re = 0.0;
      y[i].im = 0.0;
    }
  } else {
    double nt_re;
    int N2blue;
    int i;
    int rt;
    boolean_T useRadix2;
    useRadix2 = ((static_cast<int>(varargin_1) > 0) &&
                 ((static_cast<int>(varargin_1) &
                   (static_cast<int>(varargin_1) - 1)) == 0));
    N2blue = internal::fft::FFTImplementationCallback::get_algo_sizes(
        static_cast<int>(varargin_1), useRadix2, nd2);
    nt_re = 6.2831853071795862 / static_cast<double>(nd2);
    rt = static_cast<int>(static_cast<unsigned int>(nd2) >> 1) >> 1;
    costab1q.set_size(1, rt + 1);
    costab1q[0] = 1.0;
    nd2 = rt / 2 - 1;
    for (int k{0}; k <= nd2; k++) {
      costab1q[k + 1] = std::cos(nt_re * (static_cast<double>(k) + 1.0));
    }
    i = nd2 + 2;
    for (int k{i}; k < rt; k++) {
      costab1q[k] = std::sin(nt_re * static_cast<double>(rt - k));
    }
    costab1q[rt] = 0.0;
    if (!useRadix2) {
      rt = costab1q.size(1) - 1;
      nd2 = (costab1q.size(1) - 1) << 1;
      costab.set_size(1, nd2 + 1);
      sintab.set_size(1, nd2 + 1);
      costab[0] = 1.0;
      sintab[0] = 0.0;
      sintabinv.set_size(1, nd2 + 1);
      for (int k{0}; k < rt; k++) {
        sintabinv[k + 1] = costab1q[(rt - k) - 1];
      }
      i = costab1q.size(1);
      for (int k{i}; k <= nd2; k++) {
        sintabinv[k] = costab1q[k - rt];
      }
      for (int k{0}; k < rt; k++) {
        costab[k + 1] = costab1q[k + 1];
        sintab[k + 1] = -costab1q[(rt - k) - 1];
      }
      for (int k{i}; k <= nd2; k++) {
        costab[k] = -costab1q[nd2 - k];
        sintab[k] = -costab1q[k - rt];
      }
    } else {
      rt = costab1q.size(1) - 1;
      nd2 = (costab1q.size(1) - 1) << 1;
      costab.set_size(1, nd2 + 1);
      sintab.set_size(1, nd2 + 1);
      costab[0] = 1.0;
      sintab[0] = 0.0;
      for (int k{0}; k < rt; k++) {
        costab[k + 1] = costab1q[k + 1];
        sintab[k + 1] = costab1q[(rt - k) - 1];
      }
      i = costab1q.size(1);
      for (int k{i}; k <= nd2; k++) {
        costab[k] = -costab1q[nd2 - k];
        sintab[k] = costab1q[k - rt];
      }
      sintabinv.set_size(1, 0);
    }
    if (useRadix2) {
      b_x_data.set((creal_T *)&x_data[0], x_size);
      internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
          b_x_data, static_cast<int>(varargin_1), costab, sintab, y);
      if (y.size(0) > 1) {
        nt_re = 1.0 / static_cast<double>(y.size(0));
        nd2 = y.size(0);
        for (i = 0; i < nd2; i++) {
          y[i].re = nt_re * y[i].re;
          y[i].im = nt_re * y[i].im;
        }
      }
    } else {
      double b_re_tmp;
      double nt_im;
      double re_tmp;
      int nInt2;
      int nInt2m1_tmp;
      nInt2m1_tmp =
          (static_cast<int>(varargin_1) + static_cast<int>(varargin_1)) - 1;
      wwc.set_size(nInt2m1_tmp);
      rt = 0;
      wwc[static_cast<int>(varargin_1) - 1].re = 1.0;
      wwc[static_cast<int>(varargin_1) - 1].im = 0.0;
      nInt2 = static_cast<int>(varargin_1) << 1;
      for (int k{0}; k <= nfft_tmp_tmp - 2; k++) {
        nd2 = ((k + 1) << 1) - 1;
        if (nInt2 - rt <= nd2) {
          rt += nd2 - nInt2;
        } else {
          rt += nd2;
        }
        nt_im = 3.1415926535897931 * static_cast<double>(rt) /
                static_cast<double>(static_cast<int>(varargin_1));
        i = (static_cast<int>(varargin_1) - k) - 2;
        wwc[i].re = std::cos(nt_im);
        wwc[i].im = -std::sin(nt_im);
      }
      i = nInt2m1_tmp - 1;
      for (int k{i}; k >= nfft_tmp_tmp; k--) {
        wwc[k] = wwc[(nInt2m1_tmp - k) - 1];
      }
      y.set_size(nfft_tmp_tmp);
      if (static_cast<int>(varargin_1) > x_size) {
        y.set_size(nfft_tmp_tmp);
        for (i = 0; i < nfft_tmp_tmp; i++) {
          y[i].re = 0.0;
          y[i].im = 0.0;
        }
      }
      nd2 = static_cast<int>(varargin_1);
      if (nd2 > x_size) {
        nd2 = x_size;
      }
      i = static_cast<unsigned short>(nd2);
      for (int k{0}; k < i; k++) {
        rt = (static_cast<int>(varargin_1) + k) - 1;
        nt_re = wwc[rt].re;
        nt_im = wwc[rt].im;
        re_tmp = x_data[k].im;
        b_re_tmp = x_data[k].re;
        y[k].re = nt_re * b_re_tmp + nt_im * re_tmp;
        y[k].im = nt_re * re_tmp - nt_im * b_re_tmp;
      }
      i = nd2 + 1;
      for (int k{i}; k <= nfft_tmp_tmp; k++) {
        y[k - 1].re = 0.0;
        y[k - 1].im = 0.0;
      }
      internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
          y, N2blue, costab, sintab, fv);
      internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
          wwc, N2blue, costab, sintab, b_fv);
      nd2 = fv.size(0);
      b_fv.set_size(fv.size(0));
      for (i = 0; i < nd2; i++) {
        nt_re = fv[i].re;
        nt_im = b_fv[i].im;
        re_tmp = fv[i].im;
        b_re_tmp = b_fv[i].re;
        b_fv[i].re = nt_re * b_re_tmp - re_tmp * nt_im;
        b_fv[i].im = nt_re * nt_im + re_tmp * b_re_tmp;
      }
      internal::fft::FFTImplementationCallback::r2br_r2dit_trig_impl(
          b_fv, N2blue, costab, sintabinv, fv);
      if (fv.size(0) > 1) {
        nt_re = 1.0 / static_cast<double>(fv.size(0));
        nd2 = fv.size(0);
        for (i = 0; i < nd2; i++) {
          fv[i].re = nt_re * fv[i].re;
          fv[i].im = nt_re * fv[i].im;
        }
      }
      for (int k{nfft_tmp_tmp}; k <= nInt2m1_tmp; k++) {
        double ar;
        nt_re = wwc[k - 1].re;
        nt_im = fv[k - 1].im;
        re_tmp = wwc[k - 1].im;
        b_re_tmp = fv[k - 1].re;
        ar = nt_re * b_re_tmp + re_tmp * nt_im;
        nt_re = nt_re * nt_im - re_tmp * b_re_tmp;
        if (nt_re == 0.0) {
          i = k - static_cast<int>(varargin_1);
          y[i].re = ar / static_cast<double>(static_cast<int>(varargin_1));
          y[i].im = 0.0;
        } else if (ar == 0.0) {
          i = k - static_cast<int>(varargin_1);
          y[i].re = 0.0;
          y[i].im = nt_re / static_cast<double>(static_cast<int>(varargin_1));
        } else {
          i = k - static_cast<int>(varargin_1);
          y[i].re = ar / static_cast<double>(static_cast<int>(varargin_1));
          y[i].im = nt_re / static_cast<double>(static_cast<int>(varargin_1));
        }
      }
    }
  }
}

} // namespace coder

// End of code generation (ifft.cpp)
