#include <string.h>
#include "fpga_constants.h"

extern "C" {
  void tGswFFTExternMulToTLwe(TLweSample_FPGA *_accum, const TGswSampleFFT_FPGA *_gsw) {
    #pragma HLS INTERFACE m_axi port=_accum bundle=b_accum
    #pragma HLS INTERFACE m_axi port=_gsw bundle=b_gsw

    TLweSample_FPGA accum;
    accum.current_variance = _accum->current_variance;
    tGswFFTExternMulToTLwe_load_1: for(int i=0; i<=param_k; i++) {
      tGswFFTExternMulToTLwe_load_2: for(int j=0; j<param_N; j++) {
        accum.a[i].coefsT[j] = _accum->a[i].coefsT[j];
      }
    }

    TGswSampleFFT_FPGA gsw;
    gsw.k = _gsw->k;
    gsw.l = _gsw->l;
    tGswFFTExternMulToTLwe_load_3: for(int i=0; i<(param_k + 1) * param_l; i++) {
      gsw.all_samples[i].current_variance = _gsw->all_samples[i].current_variance;
      tGswFFTExternMulToTLwe_load_4: for(int j=0; j<=param_k; j++) {
        tGswFFTExternMulToTLwe_load_5: for(int k=0; k<param_Ns2; k++) {
          gsw.all_samples[i].a[j].coefsC[k] = _gsw->all_samples[i].a[j].coefsC[k];
        }
      }
    }

    IntPolynomial deca[param_kpl];
    LagrangeHalfCPolynomial decaFFT[param_kpl];
    TLweSampleFFT_FPGA tmpa[param_kpl];
    TLweSampleFFT_FPGA tmpa_final;

    loop_decomp: for(int i=0; i<=param_k; i++) {
      #pragma HLS unroll
      tGswTorus32PolynomialDecompH(&deca[i * param_l], &accum.a[i]);
    }
    loop_ifft: for(int p=0; p<param_kpl; p++) {
      #pragma HLS unroll
      IntPolynomial_ifft(&decaFFT[p], &deca[p]);
    }
    loop_clear: for(int p=0; p<param_kpl; p++) {
      #pragma HLS unroll
      tLweFFTClear(&tmpa[p]);
    }
    loop_rto: for(int p=0; p<param_kpl; p++) {
      #pragma HLS unroll
      tLweFFTAddMulRTo(&tmpa[p], &decaFFT[p], &gsw.all_samples[p]);
    }

    for(int p=0; p<param_kpl; p++) {
      for(int k=0; k<=param_k; k++) {
        for(int i=0; i<param_Ns2; i++) {
          #pragma HLS loop_flatten
          tmpa_final.a[k].coefsC[i] += tmpa[p].a[k].coefsC[i];
        }
      }
    }

    tLweFromFFTConvert(&accum, &tmpa_final);

    _accum->current_variance = accum.current_variance;
    tGswFFTExternMulToTLwe_store_1: for(int i=0; i<=param_k; i++) {
      tGswFFTExternMulToTLwe_store_2: for(int j=0; j<param_N; j++) {
        _accum->a[i].coefsT[j] = accum.a[i].coefsT[j];
      }
    }
  }
}