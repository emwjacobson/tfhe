#include "fpga_constants.h"

extern "C" {
  void TorusPolynomial_ifft(LagrangeHalfCPolynomial* result, const TorusPolynomial* p) {
    cplx* res = result->coefsC;
    const Torus32* a = p->coefsT;

    double real_inout[_2N];
    double imag_inout[_2N];

    static const double _2pm33 = 1./double(INT64_C(1)<<33);
    int32_t* aa = (int32_t*) a;
    for (int32_t i=0; i<N; i++) real_inout[i]=aa[i]*_2pm33;
    for (int32_t i=0; i<N; i++) real_inout[N+i]=-real_inout[i];
    for (int32_t i=0; i<_2N; i++) imag_inout[i]=0;

    fft_transform_reverse(real_inout, imag_inout);

    for (int32_t i=0; i<Ns2; i++) res[i]=cplx(real_inout[2*i+1],imag_inout[2*i+1]);
  }
}