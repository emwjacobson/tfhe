#include "fpga_constants.h"

extern "C" {
  void torusPolynomialMulByXaiMinusOne(TorusPolynomial *result, int32_t a, const TorusPolynomial *source) {
    Torus32 *out = result->coefsT;
    const Torus32 *in = source->coefsT;

    if (a < param_N) {
      torusPolynomialMulByXaiMinusOne_loop_1: for(int i=0; i<param_N; i++) {
        if (i<a) {
          out[i] = -in[i - a + param_N] - in[i];
        } else {
          out[i] = in[i - a] - in[i];
        }
      }
    } else {
      const int32_t aa = a - param_N;
      torusPolynomialMulByXaiMinusOne_loop_2: for(int i=0; i<param_N; i++) {
        if (i < aa) {
          out[i] = in[i - aa + param_N] - in[i];
        } else {
          out[i] = -in[i - aa] - in[i];
        }
      }
    }
  }

  void tLweMulByXaiMinusOne(TLweSample_FPGA *result, int32_t ai, const TLweSample_FPGA *bk) {
    tLweMulByXaiMinusOne_loop_1: for (int32_t i = 0; i <= param_k; i++)
      torusPolynomialMulByXaiMinusOne(&result->a[i], ai, &bk->a[i]);
  }

  void tLweAddTo(TLweSample_FPGA *result_out, const TLweSample_FPGA *result_in, const TLweSample_FPGA *sample) {
    tLweAddTo_loop_1: for (int32_t i = 0; i <= param_k; i++) {
      tLweAddTo_loop_2: for (int32_t j = 0; j < param_N; j++)
        result_out->a[i].coefsT[j] = result_in->a[i].coefsT[j] + sample->a[i].coefsT[j];
    }
    result_out->current_variance = result_in->current_variance + sample->current_variance;
  }

  void tfhe_MuxRotate_FFT(TLweSample_FPGA *result, const TLweSample_FPGA *accum, const TGswSampleFFT_FPGA *bki, const int32_t barai) {
    TLweSample_FPGA temp;
    // ACC = BKi*[(X^barai-1)*ACC]+ACC
    // temp = (X^barai-1)*ACC
    tLweMulByXaiMinusOne(&temp, barai, accum);
    // temp *= BKi
    tGswFFTExternMulToTLwe(&temp, bki);
    // ACC += temp
    tLweAddTo(result, &temp, accum);
  }
}