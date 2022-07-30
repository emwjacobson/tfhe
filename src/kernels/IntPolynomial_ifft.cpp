#include "fpga_constants.h"

extern "C" {
  void IntPolynomial_ifft(LagrangeHalfCPolynomial_Collapsed result, const IntPolynomial_Collapsed p) {
    double real_inout[_2N];
    double imag_inout[_2N];

    for (int32_t i=0; i<N; i++) real_inout[i]=p[i]/2.;
    for (int32_t i=0; i<N; i++) real_inout[N+i]=-real_inout[i];
    for (int32_t i=0; i<_2N; i++) imag_inout[i]=0;

    fft_transform_reverse(real_inout, imag_inout);

    double* res_dbl = (double*)result;
    for (int32_t i=0; i<N; i+=2) {
        res_dbl[i]=real_inout[i+1];
        res_dbl[i+1]=imag_inout[i+1];
    }
  }
}