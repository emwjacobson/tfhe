#include "fpga_constants.h"

extern "C" {
    void LagrangeHalfCPolynomialAddMul(
        LagrangeHalfCPolynomial *accum,
        const LagrangeHalfCPolynomial *a,
        const LagrangeHalfCPolynomial *b)
    {
    #pragma HLS inline
        const cplx* aa = a->coefsC;
        const cplx* bb = b->coefsC;
        cplx* rr = accum->coefsC;
        for (int32_t i=0; i<param_Ns2; i++)
            rr[i] += aa[i]*bb[i];
    }

    void tLweFFTAddMulRTo(TLweSampleFFT_FPGA *tmpa, const LagrangeHalfCPolynomial *decaFFT, const TGswSampleFFT_FPGA *gsw) {
    #pragma HLS inline
        for(int p=0; p<param_kpl; p++) {
            const LagrangeHalfCPolynomial *_p = &decaFFT[p];

            for (int32_t i = 0; i <= param_k; i++) {
                LagrangeHalfCPolynomialAddMul(&tmpa->a[i], _p, &gsw->all_samples[p].a[i]);
            }
        }
    }
}