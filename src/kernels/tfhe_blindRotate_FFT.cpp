#include "fpga_constants.h"

extern "C" {
  void tLweCopy(TLweSample_FPGA *result, const TLweSample_FPGA *sample) {
    for (int32_t i = 0; i <= param_k; ++i)
      for (int32_t j = 0; j < param_N; ++j)
        result->a[i].coefsT[j] = sample->a[i].coefsT[j];

    result->current_variance = sample->current_variance;
  }

  /**
   * multiply the accumulator by X^sum(bara_i.s_i)
   * @param accum the TLWE sample to multiply
   * @param bk An array of n TGSW FFT samples where bk_i encodes s_i
   * @param bara An array of n coefficients between 0 and 2N-1
   */
  void tfhe_blindRotate_FFT(TLweSample_FPGA *_accum, const TGswSampleFFT_FPGA bkFFT[param_n], const int32_t bara[param_n]) {
    #pragma HLS INTERFACE m_axi port=_accum bundle=b_accum
    #pragma HLS INTERFACE m_axi port=bkFFT bundle=b_bkFFT
    #pragma HLS INTERFACE m_axi port=bara bundl=b_bara

    // Copy data in
    TLweSample_FPGA accum;
    accum.current_variance = _accum->current_variance;
    for(int i=0; i<=param_k; i++) {
      for(int j=0; j<param_N; j++) {
        accum.a[i].coefsT[j] = _accum->a[i].coefsT[j];
      }
    }

    for (int32_t i = 0; i < param_n; i++) {
      const int32_t barai = bara[i];
      if (barai == 0) continue; //indeed, this is an easy case!
      tfhe_MuxRotate_FFT(&accum, &bkFFT[i], barai);
    }

    // Copy data out
    _accum->current_variance = accum.current_variance;
    for(int i=0; i<=param_k; i++) {
      for(int j=0; j<param_N; j++) {
        _accum->a[i].coefsT[j] = accum.a[i].coefsT[j];
      }
    }

    return;
  }
}