#include "fpga_constants.h"

extern "C" {
  void IntPolynomial_ifft_loop(LagrangeHalfCPolynomial *decaFFT, IntPolynomial *deca) {
    for(int p=0; p<param_kpl; p++) {
      #pragma HLS pipeline off
      LagrangeHalfCPolynomial *result = &decaFFT[p];
      IntPolynomial *_p = &deca[p];

      double real_inout[param_2N];
      double imag_inout[param_2N];

      #pragma HLS array_partition variable=real_inout type=complete
      #pragma HLS array_partition variable=imag_inout type=complete

      for (int32_t i=0; i<param_N; i++) real_inout[i]=_p->coefs[i]/2.;
      for (int32_t i=0; i<param_N; i++) real_inout[param_N+i]=-real_inout[i];
      for (int32_t i=0; i<param_2N; i++) imag_inout[i]=0;

      fft_transform_reverse(real_inout, imag_inout);

      // double* res_dbl = (double*)result->coefsC;
      // for (int32_t i=0; i<param_N; i+=2) { // i = 0, 2, 4, 6 ... 1020, 1022, 1024
      //   res_dbl[i]=real_inout[i+1]; // i+1 = 1, 3, 5, 7 ... 1021, 1023, 1025
      //   res_dbl[i+1]=imag_inout[i+1];
      // }

      for(int i=0; i<param_Ns2; i++) { // i = 0, 1, 2, 3, ... 510, 511, 512
        result->coefsC[i] = cplx(real_inout[2*i+1], imag_inout[2*i+1]); // 2*i + 1 = 1, 3, 5, 7 ... 1021, 1023, 1025
      }
    }
  }
}