#pragma once


// #include <complex>
#include <hls_x_complex.h>
#include <ap_fixed.h>


//typedef float T;
typedef ap_fixed<48, 24> T_48_24;
typedef ap_fixed<26, 2> T_26_2;
typedef hls::x_complex<T_48_24> Complex;
//using namespace std;
