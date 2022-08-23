#include "fpga_constants.h"

extern "C" {
    void LagrangeHalfCPolynomialAddMul_2(
        LagrangeHalfCPolynomial *accum,
        LagrangeHalfCPolynomial *a,
        LagrangeHalfCPolynomial *b)
    {
        cplx* aa = a->coefsC;
        cplx* bb = b->coefsC;
        cplx* rr = accum->coefsC;
        for (int32_t i=0; i<param_Ns2; i++)
            rr[i] += aa[i]*bb[i];
    }

    void tLweFFTAddMulRTo_loop(TLweSampleFFT_FPGA *tmpa, LagrangeHalfCPolynomial *decaFFT, TGswSampleFFT_FPGA *gsw) {
        for(int p=0; p<param_kpl; p++) {
            LagrangeHalfCPolynomial *_p = &decaFFT[p];

            for (int32_t i = 0; i <= param_k; i++) {
                LagrangeHalfCPolynomialAddMul_2(&tmpa->a[i], _p, &gsw->all_samples[p].a[i]);
            }
        }
    }
}