#include "fpga_constants.h"

extern "C" {
  void TorusPolynomial_ifft(LagrangeHalfCPolynomial *result, TorusPolynomial *p) {
    double real_inout[param_2N];
    double imag_inout[param_2N];

    static const double _2pm33 = 1./double(INT64_C(1)<<33);
    for (int32_t i=0; i<param_N; i++) real_inout[i]=p->coefsT[i]*_2pm33;
    for (int32_t i=0; i<param_N; i++) real_inout[param_N+i]=-real_inout[i];
    for (int32_t i=0; i<param_2N; i++) imag_inout[i]=0;

    fft_transform_reverse(real_inout, imag_inout);

    for (int32_t i=0; i<param_Ns2; i++) result->coefsC[i]=cplx(real_inout[2*i+1],imag_inout[2*i+1]);
  }
}