#include "fpga_constants.h"

extern "C" {
    void LagrangeHalfCPolynomialAddMul(
        LagrangeHalfCPolynomial *accum,
        const LagrangeHalfCPolynomial *a,
        const LagrangeHalfCPolynomial *b)
    {
        const cplx* aa = a->coefsC;
        const cplx* bb = b->coefsC;
        cplx* rr = accum->coefsC;
        for (int32_t i=0; i<param_Ns2; i++)
            rr[i] += aa[i]*bb[i];
    }

    void tLweFFTAddMulRTo(TLweSampleFFT_FPGA *result, const LagrangeHalfCPolynomial *p, const TLweSampleFFT_FPGA *sample) {
        for (int32_t i = 0; i <= param_k; i++)
            LagrangeHalfCPolynomialAddMul(&result->a[i], p, &sample->a[i]);
    }
}