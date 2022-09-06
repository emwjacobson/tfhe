#include "fpga_constants.h"

extern "C" {
  void LagrangeHalfCPolynomialClear(LagrangeHalfCPolynomial *reps) {
    LagrangeHalfCPolynomialClear_loop_1: for (int32_t i=0; i<param_Ns2; i++) {
      reps->coefsC[i] = cplx(0, 0);
    }
  }

  void tLweFFTClear_dataflow(TLweSampleFFT_FPGA *result) {
    tLweFFTClear_loop_1: for (int32_t i = 0; i <= param_k; ++i) {
      LagrangeHalfCPolynomialClear(&result->a[i]);
    }
    result->current_variance = 0.0;
  }

  void tLweFFTClear(TLweSampleFFT_FPGA result[param_kpl]) {
    #pragma HLS dataflow
    // for(int p=0; p<param_kpl; p++) {
    tLweFFTClear_dataflow(&result[0]);
    tLweFFTClear_dataflow(&result[1]);
    tLweFFTClear_dataflow(&result[2]);
    tLweFFTClear_dataflow(&result[3]);
    tLweFFTClear_dataflow(&result[4]);
    tLweFFTClear_dataflow(&result[5]);
    // }
  }
}