#include "fpga_constants.h"

extern "C" {
  void TorusPolynomial_ifft(LagrangeHalfCPolynomial *result, const TorusPolynomial *p) {
    double real_inout[param_2N];
    double imag_inout[param_2N];

    // #pragma HLS array_partition variable=real_inout cyclic factor=128
    // #pragma HLS array_partition variable=imag_inout cyclic factor=128

    static const double _2pm33 = 1./double(INT64_C(1)<<33);

    TorusPolynomial_ifft_loop_1: for (int32_t i=0; i<param_N; i++) {
      real_inout[i]=p->coefsT[i]*_2pm33;
      real_inout[param_N+i]=-real_inout[i];
      imag_inout[i]=0;
      imag_inout[param_N+i]=0;
    }

    fft_transform_reverse(real_inout, imag_inout);

    TorusPolynomial_ifft_loop_2: for (int32_t i=0; i<param_Ns2; i++) {
      result->coefsC[i]=cplx(real_inout[2*i+1],imag_inout[2*i+1]);
    }
  }
}