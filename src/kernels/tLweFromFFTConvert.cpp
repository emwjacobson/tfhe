#include "fpga_constants.h"

extern "C" {
  void TorusPolynomial_fft_dataflow(TorusPolynomial *tp, const LagrangeHalfCPolynomial *lp) {
    #pragma HLS dataflow
    // for(int i=0; i<param_k; i++) {
    TorusPolynomial_fft(&tp[0], &lp[0]);
    TorusPolynomial_fft(&tp[1], &lp[1]);
    // }
  }

  void tLweFromFFTConvert(TLweSample_FPGA *result, const TLweSampleFFT_FPGA *source) {
    TorusPolynomial_fft_dataflow(result->a, source->a);
    result->current_variance = source->current_variance;
  }
}