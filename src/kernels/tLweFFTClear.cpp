#include "fpga_constants.h"

extern "C" {
    void LagrangeHalfCPolynomialClear(LagrangeHalfCPolynomial *reps) {
        for (int32_t i=0; i<param_Ns2; i++) {
            #pragma HLS pipeline off
	        reps->coefsC[i] = cplx(0, 0);
        }
    }

    void tLweFFTClear(TLweSampleFFT_FPGA *result) {
        for (int32_t i = 0; i <= param_k; ++i) {
            #pragma HLS pipeline off
            LagrangeHalfCPolynomialClear(&result->a[i]);
        }
        result->current_variance = 0.0;
    }
}