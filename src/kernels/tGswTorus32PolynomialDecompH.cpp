#include "fpga_constants.h"

extern "C" {
  void tGswTorus32PolynomialDecompH(IntPolynomial_Collapsed *result, const TorusPolynomial_Collapsed sample) {
    uint32_t *buf = (uint32_t *)sample;

    //First, add offset to everyone
    for (int32_t j = 0; j < N; ++j) buf[j] += offset;

    //then, do the decomposition (in parallel)
    for (int32_t p = 0; p < param_l; ++p) {
        const int32_t decal = (32 - (p + 1) * Bgbit);
        int32_t *res_p = result[p];
        for (int32_t j = 0; j < N; ++j) {
            uint32_t temp1 = (buf[j] >> decal) & maskMod;
            res_p[j] = temp1 - halfBg;
        }
    }

    //finally, remove offset to everyone
    for (int32_t j = 0; j < N; ++j) buf[j] -= offset;
  }
}