#include "fpga_constants.h"

extern "C" {
  void TorusPolynomial_ifft(LagrangeHalfCPolynomial_Collapsed result, const TorusPolynomial_Collapsed p) {
    double real_inout[_2N];
    double imag_inout[_2N];

    static const double _2pm33 = 1./double(INT64_C(1)<<33);
    for (int32_t i=0; i<N; i++) real_inout[i]=p[i]*_2pm33;
    for (int32_t i=0; i<N; i++) real_inout[N+i]=-real_inout[i];
    for (int32_t i=0; i<_2N; i++) imag_inout[i]=0;

    fft_transform_reverse(real_inout, imag_inout);

    for (int32_t i=0; i<Ns2; i++) result[i]=cplx(real_inout[2*i+1],imag_inout[2*i+1]);
  }
}