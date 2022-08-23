#include "fpga_constants.h"

extern "C" {
    void tLweFromFFTConvert(TLweSample_FPGA *result, const TLweSampleFFT_FPGA *source) {
        for (int32_t i = 0; i <= param_k; ++i) {
            TorusPolynomial_fft(&result->a[i], &source->a[i]);
        }
        result->current_variance = source->current_variance;
    }
}