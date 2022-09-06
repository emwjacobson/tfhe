#include "fpga_constants.h"

extern "C" {
  void IntPolynomial_ifft_dataflow(LagrangeHalfCPolynomial *result, const IntPolynomial *p) {
    double real_inout[param_2N];
    double imag_inout[param_2N];

    // #pragma HLS array_partition variable=real_inout cyclic factor=512
    // #pragma HLS array_partition variable=imag_inout cyclic factor=512
    // // #pragma HLS bind_storage variable=real_inout type=RAM_T2P
    // // #pragma HLS bind_storage variable=imag_inout type=RAM_T2P

    IntPolynomial_ifft_loop_1: for (int32_t i=0; i<param_N; i++) {
      real_inout[i] = p->coefs[i]/2.;
      real_inout[param_N+i] = -real_inout[i];
    }

    IntPolynomial_ifft_loop_2: for(int32_t i=0; i<param_2N; i++) {
      imag_inout[i] = 0;
    }

    fft_transform(imag_inout, real_inout);

    IntPolynomial_ifft_loop_3: for(int i=0; i<param_Ns2; i++) { // i = 0, 1, 2, 3, ... 510, 511, 512
      result->coefsC[i] = cplx(real_inout[2*i+1], imag_inout[2*i+1]); // 2*i + 1 = 1, 3, 5, 7 ... 1021, 1023, 1025
    }
  }

  void IntPolynomial_ifft(LagrangeHalfCPolynomial result[param_kpl], const IntPolynomial p[param_kpl]) {
    #pragma HLS dataflow
    // for(int p=0; p<param_kpl; p++) {
    IntPolynomial_ifft_dataflow(&result[0], &p[0]);
    IntPolynomial_ifft_dataflow(&result[1], &p[1]);
    IntPolynomial_ifft_dataflow(&result[2], &p[2]);
    IntPolynomial_ifft_dataflow(&result[3], &p[3]);
    IntPolynomial_ifft_dataflow(&result[4], &p[4]);
    IntPolynomial_ifft_dataflow(&result[5], &p[5]);
    // }
  }
}