#include "fpga_constants.h"

extern "C" {
  void TorusPolynomial_fft(TorusPolynomial *result, const LagrangeHalfCPolynomial *p) {
    double real_inout[param_2N];
    double imag_inout[param_2N];

    // #pragma HLS array_partition variable=real_inout cyclic factor=128
    // #pragma HLS array_partition variable=imag_inout cyclic factor=128

    //double* a_dbl=(double*) a;
    TorusPolynomial_fft_loop_1: for (int32_t i=0; i<param_N; i++){
      real_inout[2*i] = 0;
      imag_inout[2*i] = 0;
    }
    TorusPolynomial_fft_loop_2: for (int32_t i=0; i<param_Ns2; i++) {
      real_inout[2*i+1] = real(p->coefsC[i]);
      imag_inout[2*i+1] = imag(p->coefsC[i]);
      real_inout[param_2N-1-2*i] = real(p->coefsC[i]);
      imag_inout[param_2N-1-2*i] = -imag(p->coefsC[i]);
    }
    fft_transform(real_inout, imag_inout);

    static const double _2p32 = double(INT64_C(1)<<32);
    static const double _1sN = double(1)/double(param_N);

    TorusPolynomial_fft_loop_3: for (int32_t i=0; i<param_N; i++) {
      result->coefsT[i] = (Torus32)((int64_t)(real_inout[i] * _1sN * _2p32));
    }
  }
}