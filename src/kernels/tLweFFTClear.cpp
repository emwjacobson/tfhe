#include "fpga_constants.h"

extern "C" {
    void LagrangeHalfCPolynomialClear(LagrangeHalfCPolynomial *reps) {
    #pragma HLS inline
        for (int32_t i=0; i<param_Ns2; i++) {
	        reps->coefsC[i] = cplx(0, 0);
        }
    }

    void tLweFFTClear(TLweSampleFFT_FPGA *result) {
    #pragma HLS inline
        for (int32_t i = 0; i <= param_k; ++i) {
            LagrangeHalfCPolynomialClear(&result->a[i]);
        }
        result->current_variance = 0.0;
    }
}