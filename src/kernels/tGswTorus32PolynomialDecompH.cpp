#include "fpga_constants.h"

extern "C" {
  void tGswTorus32PolynomialDecompH_dataflow(IntPolynomial result[param_l], const TorusPolynomial *sample) {
    uint32_t buf[param_N];

    //First, add offset to everyone
    tGswTorus32PolynomialDecompH_loop_1: for (int32_t j = 0; j < param_N; ++j) {
      buf[j] = (uint32_t)(sample->coefsT[j] + param_offset);
    }

    //then, do the decomposition (in parallel)
    tGswTorus32PolynomialDecompH_loop_2: for (int32_t p = 0; p < param_l; ++p){
      const int32_t decal = (32 - ((p + 1) * param_Bgbit));
      tGswTorus32PolynomialDecompH_loop_3: for (int32_t j = 0; j < param_N; ++j){
        uint32_t temp1 = (buf[j] >> decal) & param_maskMod;
        result[p].coefs[j] = temp1 - param_halfBg;
      }
    }
  }

  void tGswTorus32PolynomialDecompH(IntPolynomial result[param_kpl], const TLweSample_FPGA *accum) {
    #pragma HLS dataflow
    // for(int i=0; i<=param_k; i++) {
    tGswTorus32PolynomialDecompH_dataflow(&result[0], &accum->a[0]);
    tGswTorus32PolynomialDecompH_dataflow(&result[3], &accum->a[1]);
    // }
  }
}