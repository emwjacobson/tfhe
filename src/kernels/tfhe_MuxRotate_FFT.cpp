#include "fpga_constants.h"

extern "C" {
  void torusPolynomialMulByXaiMinusOne(TorusPolynomial *result, int32_t a, const TorusPolynomial *source) {
    Torus32 *out = result->coefsT;
    const Torus32 *in = source->coefsT;

    if (a < param_N) {
      for (int32_t i = 0; i < a; i++)//sur que i-a<0
        out[i] = -in[i - a + param_N] - in[i];
      for (int32_t i = a; i < param_N; i++)//sur que N>i-a>=0
        out[i] = in[i - a] - in[i];
    } else {
      const int32_t aa = a - param_N;
      for (int32_t i = 0; i < aa; i++)//sur que i-a<0
        out[i] = in[i - aa + param_N] - in[i];
      for (int32_t i = aa; i < param_N; i++)//sur que N>i-a>=0
        out[i] = -in[i - aa] - in[i];
    }
  }

  void tLweMulByXaiMinusOne(TLweSample_FPGA *result, int32_t ai, const TLweSample_FPGA *bk) {
    for (int32_t i = 0; i <= param_k; i++)
      torusPolynomialMulByXaiMinusOne(&result->a[i], ai, &bk->a[i]);
  }

  void torusPolynomialAddTo(TorusPolynomial *result, const TorusPolynomial *poly2) {
    Torus32 *r = result->coefsT;
    const Torus32 *b = poly2->coefsT;

    for (int32_t i = 0; i < param_N; ++i)
      r[i] += b[i];
  }

  void tLweAddTo(TLweSample_FPGA *result, const TLweSample_FPGA *sample) {
    for (int32_t i = 0; i < param_k; ++i)
      torusPolynomialAddTo(&result->a[i], &sample->a[i]);
    torusPolynomialAddTo(&result->a[param_k], &sample->a[param_k]);
    result->current_variance += sample->current_variance;
  }

  void tfhe_MuxRotate_FFT(TLweSample_FPGA *accum, const TGswSampleFFT_FPGA *bki, const int32_t barai) {
    TLweSample_FPGA temp;

    // ACC = BKi*[(X^barai-1)*ACC]+ACC
    // temp = (X^barai-1)*ACC
    tLweMulByXaiMinusOne(&temp, barai, accum);
    // temp *= BKi
    tGswFFTExternMulToTLwe(&temp, bki);
    // ACC += temp
    tLweAddTo(accum, &temp);
  }
}