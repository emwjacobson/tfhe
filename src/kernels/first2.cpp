#include <string.h>
#include "fpga_constants.h"

extern "C" {
  void first2(LagrangeHalfCPolynomial *_decaFFT, IntPolynomial *_deca, TLweSample_FPGA *_accum) {
    #pragma HLS INTERFACE m_axi port=_decaFFT bundle=b_decaFFT
    #pragma HLS INTERFACE m_axi port=_accum bundle=b_accum

    TLweSample_FPGA accum;
    accum.current_variance = _accum->current_variance;
    for(int i=0; i<=param_k; i++) {
      for(int j=0; j<param_N; j++) {
        accum.a[i].coefsT[j] = _accum->a[i].coefsT[j];
      }
    }

    IntPolynomial deca[param_kpl];
    LagrangeHalfCPolynomial decaFFT[param_kpl];

    tGswTorus32PolynomialDecompH_l(deca, &accum);
    IntPolynomial_ifft_loop(decaFFT, deca);

    for(int i=0; i<param_kpl; i++) {
      for(int j=0; j<param_N; j++) {
        _deca[i].coefs[j] = deca[i].coefs[j];
      }
    }

    for(int i=0; i<param_kpl; i++) {
      for(int j=0; j<param_Ns2; j++) {
        _decaFFT[i].coefsC[j] = decaFFT[i].coefsC[j];
      }
    }
  }
}