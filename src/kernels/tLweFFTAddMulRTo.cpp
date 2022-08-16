#include "fpga_constants.h"

extern "C" {
    void LagrangeHalfCPolynomialAddMul(
        LagrangeHalfCPolynomial_Collapsed accum,
        const LagrangeHalfCPolynomial_Collapsed a,
        const LagrangeHalfCPolynomial_Collapsed b)
    {
        const cplx* aa = a;
        const cplx* bb = b;
        cplx* rr = accum;
        for (int32_t i=0; i<param_Ns2; i++)
            rr[i] += aa[i]*bb[i];
    }

    void tLweFFTAddMulRTo(TLweSampleFFT_FPGA *result, const LagrangeHalfCPolynomial_Collapsed p, const TLweSampleFFT_FPGA *sample) {
        for (int32_t i = 0; i <= param_k; i++)
            LagrangeHalfCPolynomialAddMul(result->a[i], p, sample->a[i]);
    }
}