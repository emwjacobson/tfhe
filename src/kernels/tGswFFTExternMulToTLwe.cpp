#include <string.h>
#include "fpga_constants.h"

extern "C" {
  void tGswFFTExternMulToTLwe(TLweSample_FPGA *_accum, const TGswSampleFFT_FPGA *_gsw) {
    #pragma HLS INTERFACE m_axi port=_accum bundle=b_accum
    #pragma HLS INTERFACE m_axi port=_gsw bundle=b_gsw

    IntPolynomial_Collapsed deca[param_kpl];
    LagrangeHalfCPolynomial_Collapsed decaFFT[param_kpl];
    TLweSampleFFT_FPGA tmpa;

    for(int i=0; i<=param_k; i++) {
      tGswTorus32PolynomialDecompH(deca + i * param_l, _accum->a[i]);
    }
    for(int p=0; p<param_kpl; p++) {
      IntPolynomial_ifft(decaFFT[p], deca[p]);
    }
    tLweFFTClear(&tmpa);
    for(int p=0; p<param_kpl; p++) {
      tLweFFTAddMulRTo(&tmpa, decaFFT[p], _gsw->all_samples + p);
    }
    tLweFromFFTConvert(_accum, &tmpa);






    // #pragma HLS INTERFACE m_axi port=_accum bundle=b_accum
    // #pragma HLS INTERFACE m_axi port=_gsw bundle=b_gsw

    // // Copy data in
    // TLweSample_FPGA accum;
    // accum.current_variance = _accum->current_variance;
    // for (int i = 0; i <= param_k; i++) {
    //   for (int j = 0; j < param_N; j++) {
    //     accum.a[i][j] = _accum->a[i][j];
    //   }
    // }

    // TGswSampleFFT_FPGA gsw;
    // gsw.k = _gsw->k;
    // gsw.l = _gsw->l;
    // for (int i = 0; i < (param_k+1) * param_l; i++) {
    //   gsw.all_samples[i].current_variance = _gsw->all_samples[i].current_variance;
    //   for (int j = 0; j <= param_k; j++){
    //     for (int k = 0; k < param_Ns2; k++) {
    //       gsw.all_samples[i].a[j][k] = _gsw->all_samples[i].a[j][k];
    //     }
    //   }
    // }

    // // Local Data
    // IntPolynomial_Collapsed deca[param_kpl];
    // LagrangeHalfCPolynomial_Collapsed decaFFT[param_kpl];
    // TLweSampleFFT_FPGA tmpa;

    // // Calculations
    // for (int32_t i = 0; i <= param_k; i++)
    //   tGswTorus32PolynomialDecompH(deca + i * param_l, accum.a[i]);

    // for (int32_t p = 0; p < param_kpl; p++)
    //   IntPolynomial_ifft(decaFFT[p], deca[p]);

    // tLweFFTClear(&tmpa);
    // for (int32_t p = 0; p < param_kpl; p++) {
    //   tLweFFTAddMulRTo(&tmpa, decaFFT[p], &gsw.all_samples[p]);
    // }
    // tLweFromFFTConvert(&accum, &tmpa);

    // // Copy data out
    // _accum->current_variance = accum.current_variance;
    // for (int i = 0; i <= param_k; i++) {
    //   for (int j = 0; j < param_N; j++) {
    //     _accum->a[i][j] = accum.a[i][j];
    //   }
    // }
  }
}