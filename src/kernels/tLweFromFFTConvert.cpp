#include "fpga_constants.h"

extern "C" {
  void tLweFromFFTConvert(TLweSample_FPGA *result, const TLweSampleFFT_FPGA *source) {
    tLweFromFFTConvert_loop_1: for (int32_t i = 0; i <= param_k; ++i) {
      #pragma HLS unroll
      TorusPolynomial_fft(&result->a[i], &source->a[i]);
    }
    result->current_variance = source->current_variance;
  }
}