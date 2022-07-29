#include "fpga_constants.h"

extern "C" {
  void IntPolynomial_ifft(LagrangeHalfCPolynomial* result, const IntPolynomial* p) {
    const int32_t* a = p->coefs;

    double real_inout[_2N];
    double imag_inout[_2N];

    for (int32_t i=0; i<N; i++) real_inout[i]=a[i]/2.;
    for (int32_t i=0; i<N; i++) real_inout[N+i]=-real_inout[i];
    for (int32_t i=0; i<_2N; i++) imag_inout[i]=0;

    fft_transform_reverse(real_inout, imag_inout);

    double* res_dbl=(double*) result->coefsC;
    for (int32_t i=0; i<N; i+=2) {
        res_dbl[i]=real_inout[i+1];
        res_dbl[i+1]=imag_inout[i+1];
    }
  }
}