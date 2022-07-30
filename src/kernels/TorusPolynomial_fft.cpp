#include "fpga_constants.h"

extern "C" {
  void TorusPolynomial_fft(TorusPolynomial_Collapsed result, const LagrangeHalfCPolynomial_Collapsed p) {
    double real_inout[_2N];
    double imag_inout[_2N];

    static const double _2p32 = double(INT64_C(1)<<32);
    static const double _1sN = double(1)/double(N);
    //double* a_dbl=(double*) a;
    for (int32_t i=0; i<N; i++) real_inout[2*i]=0;
    for (int32_t i=0; i<N; i++) imag_inout[2*i]=0;
    for (int32_t i=0; i<Ns2; i++) real_inout[2*i+1]=real(p[i]);
    for (int32_t i=0; i<Ns2; i++) imag_inout[2*i+1]=imag(p[i]);
    for (int32_t i=0; i<Ns2; i++) real_inout[_2N-1-2*i]=real(p[i]);
    for (int32_t i=0; i<Ns2; i++) imag_inout[_2N-1-2*i]=-imag(p[i]);
    fft_transform(real_inout, imag_inout);
    for (int32_t i=0; i<N; i++) result[i]=Torus32(int64_t(real_inout[i]*_1sN*_2p32));
    //pas besoin du fmod... Torus32(int64_t(fmod(rev_out[i]*_1sN,1.)*_2p32));
  }
}