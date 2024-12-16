#include "mat_utils.h"
#include "jacobi_evd.h"
#include <hls_stream.h>
#include "ap_axi_sdata.h"

typedef hls::axis<double, 0, 0, 1> axis_double;


void jacobi_svd(hls::stream<axis_double>& A_in_stream, hls::stream<axis_double>& U_out_stream, hls::stream<axis_double>& V_out_stream, hls::stream<axis_double>& S_out_stream, int m, int n);
