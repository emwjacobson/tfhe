#include <string.h>
#include "fpga_constants.h"

extern "C" {
  void tGswFFTExternMulToTLwe(TLweSample_FPGA *_accum, TGswSampleFFT_FPGA *_gsw) {
    #pragma HLS INTERFACE m_axi port=_accum bundle=b_accum
    #pragma HLS INTERFACE m_axi port=_gsw bundle=b_gsw

    TLweSample_FPGA accum;
    accum.current_variance = _accum->current_variance;
    for(int i=0; i<=param_k; i++) {
      for(int j=0; j<param_N; j++) {
        accum.a[i].coefsT[j] = _accum->a[i].coefsT[j];
      }
    }

    TGswSampleFFT_FPGA gsw;
    gsw.k = _gsw->k;
    gsw.l = _gsw->l;
    for(int i=0; i<(param_k + 1) * param_l; i++) {
      gsw.all_samples[i].current_variance = _gsw->all_samples[i].current_variance;
      for(int j=0; j<=param_k; j++) {
        for(int k=0; k<param_Ns2; k++) {
          gsw.all_samples[i].a[j].coefsC[k] = _gsw->all_samples[i].a[j].coefsC[k];
        }
      }
    }

    IntPolynomial deca[param_kpl];
    LagrangeHalfCPolynomial decaFFT[param_kpl];
    TLweSampleFFT_FPGA tmpa;

    #pragma HLS dataflow
    tGswTorus32PolynomialDecompH_l(deca, &accum);
    IntPolynomial_ifft_loop(decaFFT, deca);
    tLweFFTClear(&tmpa);
    tLweFFTAddMulRTo_loop(&tmpa, decaFFT, &gsw);
    tLweFromFFTConvert(&accum, &tmpa);

    _accum->current_variance = accum.current_variance;
    for(int i=0; i<=param_k; i++) {
      for(int j=0; j<param_N; j++) {
        _accum->a[i].coefsT[j] = accum.a[i].coefsT[j];
      }
    }
  }
}