#include "fpga_constants.h"
#include <assert.h>

extern "C" {
  void tLweCopy(TLweSample_FPGA *result, const TLweSample_FPGA *sample) {
    tLweCopy_loop_1: for (int32_t i = 0; i <= param_k; ++i)
      tLweCopy_loop_2: for (int32_t j = 0; j < param_N; ++j)
        result->a[i].coefsT[j] = sample->a[i].coefsT[j];
    result->current_variance = sample->current_variance;
  }


  void tfhe_blindRotate_FFT(TLweSample_FPGA *_accum, const TGswSampleFFT_FPGA _bkFFT[param_n], const int32_t _bara[param_n]) {
    #pragma HLS INTERFACE m_axi port=_accum bundle=b_accum
    #pragma HLS INTERFACE m_axi port=_bkFFT bundle=b_bkFFT
    #pragma HLS INTERFACE m_axi port=_bara bundle=b_bara

    TLweSample_FPGA accum;
    accum.current_variance = _accum->current_variance;
    tfhe_blindRotate_FFT_load_1: for(int i=0; i<=param_k; i++) {
      tfhe_blindRotate_FFT_load_2: for(int j=0; j<param_N; j++) {
        accum.a[i].coefsT[j] = _accum->a[i].coefsT[j];
      }
    }

    TLweSample_FPGA temp;
    TLweSample_FPGA *temp2 = &temp;
    TLweSample_FPGA *temp3 = &accum;

    for (int32_t i = 0; i < param_n; i++) {
      const int32_t barai = _bara[i];
      // if (barai == 0) continue; //indeed, this is an easy case!

      tfhe_MuxRotate_FFT(temp2, temp3, &_bkFFT[i], barai);

      // swap(temp2, temp3);
      TLweSample_FPGA *t = temp2;
      temp2 = temp3;
      temp3 = t;
    }

    assert(param_n % 2 == 0);
    // if (temp3 != &accum) {
    //   tLweCopy(accum, temp3);
    // }

    _accum->current_variance = accum.current_variance;
    tfhe_blindRotate_FFT_store_1: for(int i=0; i<=param_k; i++) {
      tfhe_blindRotate_FFT_store_2: for(int j=0; j<param_N; j++) {
        _accum->a[i].coefsT[j] = accum.a[i].coefsT[j];
      }
    }
  }
}