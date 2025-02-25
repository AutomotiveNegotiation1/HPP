//
// findpeaks.cpp
//
// Code generation for function 'findpeaks'
//

// Include files
#include "findpeaks.h"
#include "eml_setop.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>

// Function Definitions
namespace coder {
void findpeaks(const array<double, 2U> &Yin, array<double, 2U> &Ypk,
               array<double, 2U> &Xpk)
{
  array<int, 2U> y;
  array<int, 1U> c;
  array<int, 1U> fPk;
  array<int, 1U> iInfinite;
  array<int, 1U> iPk;
  array<int, 1U> idx;
  double ykfirst;
  int i;
  int kfirst;
  int nInf;
  int nPk;
  char dir;
  boolean_T isinfykfirst;
  i = Yin.size(1);
  fPk.set_size(Yin.size(1));
  iInfinite.set_size(Yin.size(1));
  nPk = 0;
  nInf = 0;
  dir = 'n';
  kfirst = 0;
  ykfirst = rtInf;
  isinfykfirst = true;
  for (int k{1}; k <= i; k++) {
    double yk;
    boolean_T isinfyk;
    yk = Yin[k - 1];
    if (std::isnan(yk)) {
      yk = rtInf;
      isinfyk = true;
    } else if (std::isinf(yk) && (yk > 0.0)) {
      isinfyk = true;
      nInf++;
      iInfinite[nInf - 1] = k;
    } else {
      isinfyk = false;
    }
    if (yk != ykfirst) {
      char previousdir;
      previousdir = dir;
      if (isinfyk || isinfykfirst) {
        dir = 'n';
      } else if (yk < ykfirst) {
        dir = 'd';
        if (previousdir == 'i') {
          nPk++;
          fPk[nPk - 1] = kfirst;
        }
      } else {
        dir = 'i';
      }
      ykfirst = yk;
      kfirst = k;
      isinfykfirst = isinfyk;
    }
  }
  if (nPk < 1) {
    i = 0;
  } else {
    i = nPk;
  }
  fPk.set_size(i);
  if (nInf < 1) {
    nInf = 0;
  }
  iInfinite.set_size(nInf);
  iPk.set_size(i);
  nPk = 0;
  for (int k{0}; k < i; k++) {
    ykfirst = Yin[fPk[k] - 1];
    if ((ykfirst > 0.9) &&
        (ykfirst - std::fmax(Yin[fPk[k] - 2], Yin[fPk[k]]) >= 0.0)) {
      nPk++;
      iPk[nPk - 1] = fPk[k];
    }
  }
  if (nPk < 1) {
    nPk = 0;
  }
  iPk.set_size(nPk);
  do_vectors(iPk, iInfinite, c, idx, fPk);
  nInf = c.size(0);
  y.set_size(1, c.size(0));
  if (c.size(0) > 0) {
    y[0] = 1;
    nPk = 1;
    for (int k{2}; k <= nInf; k++) {
      nPk++;
      y[k - 1] = nPk;
    }
  }
  idx.set_size(c.size(0));
  for (i = 0; i < nInf; i++) {
    idx[i] = y[i];
  }
  if (idx.size(0) > Yin.size(1)) {
    fPk.set_size(Yin.size(1));
    idx.set_size(Yin.size(1));
  } else {
    fPk.set_size(c.size(0));
  }
  nPk = fPk.size(0);
  fPk.set_size(nPk);
  for (i = 0; i < nPk; i++) {
    fPk[i] = c[idx[i] - 1];
  }
  Ypk.set_size(1, nPk);
  Xpk.set_size(1, nPk);
  for (i = 0; i < nPk; i++) {
    Ypk[i] = Yin[fPk[i] - 1];
    Xpk[i] = static_cast<unsigned int>(fPk[i]);
  }
}

} // namespace coder

// End of code generation (findpeaks.cpp)
