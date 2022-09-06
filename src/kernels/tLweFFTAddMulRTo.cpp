#include "fpga_constants.h"

extern "C" {
  void LagrangeHalfCPolynomialAddMul(LagrangeHalfCPolynomial *accum, const LagrangeHalfCPolynomial *a, const LagrangeHalfCPolynomial *b) {
    const cplx* aa = a->coefsC;
    const cplx* bb = b->coefsC;
    cplx* rr = accum->coefsC;
    LagrangeHalfCPolynomialAddMul_loop_1: for (int32_t i=0; i<param_Ns2; i++) {
      // accum should have been zero'd out before, and with splitting `tmpa` into multiple arrays, can avoid the +=
      // rr[i] += aa[i]*bb[i];
      rr[i] = aa[i] * bb[i];
    }
  }

  void tLweFFTAddMulRTo_dataflow(TLweSampleFFT_FPGA *result, const LagrangeHalfCPolynomial *p, const TLweSampleFFT_FPGA *sample) {
    tLweFFTAddMulRTo_loop_1: for (int32_t i = 0; i <= param_k; i++) {
      LagrangeHalfCPolynomialAddMul(&result->a[i], p, &sample->a[i]);
    }
  }

  void tLweFFTAddMulRTo(TLweSampleFFT_FPGA result[param_kpl], const LagrangeHalfCPolynomial p[param_kpl], const TGswSampleFFT_FPGA *gsw) {
    #pragma HLS dataflow
    // for(int p=0; p<param_kpl; p++) {
    tLweFFTAddMulRTo_dataflow(&result[0], &p[0], &gsw->all_samples[0]);
    tLweFFTAddMulRTo_dataflow(&result[1], &p[1], &gsw->all_samples[1]);
    tLweFFTAddMulRTo_dataflow(&result[2], &p[2], &gsw->all_samples[2]);
    tLweFFTAddMulRTo_dataflow(&result[3], &p[3], &gsw->all_samples[3]);
    tLweFFTAddMulRTo_dataflow(&result[4], &p[4], &gsw->all_samples[4]);
    tLweFFTAddMulRTo_dataflow(&result[5], &p[5], &gsw->all_samples[5]);
    // }
  }
}