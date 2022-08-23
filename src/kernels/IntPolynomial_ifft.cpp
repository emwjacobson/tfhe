#include "fpga_constants.h"

extern "C" {
  void IntPolynomial_ifft(LagrangeHalfCPolynomial *result, const IntPolynomial *p) {
    double real_inout[param_2N];
    double imag_inout[param_2N];

    for (int32_t i=0; i<param_N; i++) real_inout[i]=p->coefs[i]/2.;
    for (int32_t i=0; i<param_N; i++) real_inout[param_N+i]=-real_inout[i];
    for (int32_t i=0; i<param_2N; i++) imag_inout[i]=0;

    fft_transform_reverse(real_inout, imag_inout);

    // double* res_dbl = (double*)result->coefsC;
    // for (int32_t i=0; i<param_N; i+=2) {
    //   res_dbl[i]=real_inout[i+1];
    //   res_dbl[i+1]=imag_inout[i+1];
    // }
    for(int i=0; i<param_Ns2; i++) { // i = 0, 1, 2, 3, ... 510, 511, 512
      result->coefsC[i] = cplx(real_inout[2*i+1], imag_inout[2*i+1]); // 2*i + 1 = 1, 3, 5, 7 ... 1021, 1023, 1025
    }
  }
}