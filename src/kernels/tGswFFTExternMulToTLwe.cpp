#include <string.h>
#include "fpga_constants.h"

extern "C" {
  void tmpa_accum(TLweSampleFFT_FPGA *tmpa_final, TLweSampleFFT_FPGA tmpa[param_kpl]) {
    loop_tmpa_final_1: for(int p=0; p<param_kpl; p++) {
      loop_tmpa_final_2: for(int k=0; k<=param_k; k++) {
        loop_tmpa_final_3: for(int i=0; i<param_Ns2; i++) {
          tmpa_final->a[k].coefsC[i] += tmpa[p].a[k].coefsC[i];
        }
      }
    }
  }

  void tGswFFTExternMulToTLwe(TLweSample_FPGA *_accum, const TGswSampleFFT_FPGA *_gsw) {
    #pragma HLS INTERFACE m_axi port=_accum bundle=b_accum
    #pragma HLS INTERFACE m_axi port=_gsw bundle=b_gsw

    TLweSample_FPGA accum;
    #pragma HLS array_partition variable=accum.a complete dim=1
    accum.current_variance = _accum->current_variance;
    tGswFFTExternMulToTLwe_load_1: for(int i=0; i<=param_k; i++) {
      tGswFFTExternMulToTLwe_load_2: for(int j=0; j<param_N; j++) {
        accum.a[i].coefsT[j] = _accum->a[i].coefsT[j];
      }
    }

    TGswSampleFFT_FPGA gsw;
    #pragma HLS array_partition variable=gsw.all_samples complete dim=1
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

    #pragma HLS array_partition variable=deca complete dim=1
    #pragma HLS array_partition variable=decaFFT complete dim=1
    #pragma HLS array_partition variable=tmpa complete dim=1
    #pragma HLS array_partition variable=tmpa_final.a complete dim=1

    tGswTorus32PolynomialDecompH(deca, &accum);
    IntPolynomial_ifft(decaFFT, deca);
    tLweFFTClear(tmpa);
    tLweFFTAddMulRTo(tmpa, decaFFT, &gsw);
    tmpa_accum(&tmpa_final, tmpa);
    tLweFromFFTConvert(&accum, &tmpa_final);

    _accum->current_variance = accum.current_variance;
    tGswFFTExternMulToTLwe_store_1: for(int i=0; i<=param_k; i++) {
      tGswFFTExternMulToTLwe_store_2: for(int j=0; j<param_N; j++) {
        _accum->a[i].coefsT[j] = accum.a[i].coefsT[j];
      }
    }
  }
}