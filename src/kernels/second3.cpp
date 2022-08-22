#include <string.h>
#include "fpga_constants.h"

extern "C" {
  void second3(TLweSample_FPGA *_accum, LagrangeHalfCPolynomial *_decaFFT, TGswSampleFFT_FPGA *_gsw) {
    #pragma HLS INTERFACE m_axi port=_accum bundle=b_accum
    #pragma HLS INTERFACE m_axi port=_decaFFT bundle=b_decaFFT
    #pragma HLS INTERFACE m_axi port=_gsw bundle=b_gsw

    TLweSample_FPGA accum;
    LagrangeHalfCPolynomial decaFFT[param_kpl];
    TGswSampleFFT_FPGA gsw;

    accum.current_variance = _accum->current_variance;
    for(int i=0; i<=param_k; i++) {
      for(int j=0; j<param_N; j++) {
        accum.a[i].coefsT[j] = _accum->a[i].coefsT[j];
      }
    }

    for(int i=0; i<param_kpl; i++) {
      for(int j=0; j<param_Ns2; j++) {
        decaFFT[i].coefsC[j] = _decaFFT[i].coefsC[j];
      }
    }

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

    TLweSampleFFT_FPGA tmpa;

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