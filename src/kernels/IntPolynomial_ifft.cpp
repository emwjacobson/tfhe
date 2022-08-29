#include "fpga_constants.h"

extern "C" {
  void IntPolynomial_ifft(LagrangeHalfCPolynomial *result, const IntPolynomial *p) {
    double real_inout[param_2N];
    double imag_inout[param_2N];

    // #pragma HLS array_partition variable=real_inout cyclic factor=128
    // #pragma HLS array_partition variable=imag_inout cyclic factor=128

    IntPolynomial_ifft_loop_1: for (int32_t i=0; i<param_N; i++) {
      real_inout[i] = p->coefs[i]/2.;
      real_inout[param_N+i] = -real_inout[i];
      imag_inout[i] = 0;
      imag_inout[param_N+i] = 0;
    }

    fft_transform(imag_inout, real_inout);

    IntPolynomial_ifft_loop_2: for(int i=0; i<param_Ns2; i++) { // i = 0, 1, 2, 3, ... 510, 511, 512
      result->coefsC[i] = cplx(real_inout[2*i+1], imag_inout[2*i+1]); // 2*i + 1 = 1, 3, 5, 7 ... 1021, 1023, 1025
    }
  }
}