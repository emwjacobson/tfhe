#include "fpga_constants.h"

extern "C" {
  void TorusPolynomial_fft(TorusPolynomial *result, const LagrangeHalfCPolynomial *p) {
    double real_inout[param_2N];
    double imag_inout[param_2N];

    // #pragma HLS array_partition variable=real_inout cyclic factor=512
    // #pragma HLS array_partition variable=imag_inout cyclic factor=512
    // // #pragma HLS bind_storage variable=real_inout type=RAM_T2P
    // // #pragma HLS bind_storage variable=imag_inout type=RAM_T2P

    // Even Elements = 0
    TorusPolynomial_fft_loop_1: for (int32_t i=0; i<param_N; i++){
      real_inout[2*i] = 0;
      imag_inout[2*i] = 0;
    }
    // Odd elements
    TorusPolynomial_fft_loop_2: for (int32_t i=0; i<param_Ns2; i++) {
      real_inout[2*i+1] = real(p->coefsC[i]); // 1, 3, 5, ..., 1021, 1023
      imag_inout[2*i+1] = imag(p->coefsC[i]); // 1, 3, 5, ..., 1021, 1023
      real_inout[param_2N-1-2*i] = real(p->coefsC[i]); // 2047, 2045, 2043, ..., 1027, 1025
      imag_inout[param_2N-1-2*i] = -imag(p->coefsC[i]); // 2047, 2045, 2043, ..., 1027, 1025
    }
    fft_transform(real_inout, imag_inout);

    static const double _2p32 = double(INT64_C(1)<<32);
    static const double _1sN = double(1)/double(param_N);

    TorusPolynomial_fft_loop_3: for (int32_t i=0; i<param_N; i++) {
      result->coefsT[i] = (Torus32)((int64_t)(real_inout[i] * _1sN * _2p32));
    }
  }
}