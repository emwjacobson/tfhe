#ifndef TFHE_TEST_ENVIRONMENT
/* ***************************************************
TGSW fft operations
*************************************************** */

#include <cstdlib>
#include <iostream>
#include <random>
#include <cassert>
#include <ccomplex>
#include "tfhe_core.h"
#include "numeric_functions.h"
#include "lweparams.h"
#include "lwekey.h"
#include "lwesamples.h"
#include "lwe-functions.h"
#include "tlwe_functions.h"
#include "tgsw_functions.h"
#include "polynomials_arithmetic.h"
#include "lagrangehalfc_arithmetic.h"
#include "lwebootstrappingkey.h"

using namespace std;
#else
#undef EXPORT
#define EXPORT
#endif


//constructor content
EXPORT void init_TGswSampleFFT(TGswSampleFFT *obj, const TGswParams *params) {
    const int32_t k = params->tlwe_params->k;
    const int32_t l = params->l;
    TLweSampleFFT *all_samples = new_TLweSampleFFT_array((k + 1) * l, params->tlwe_params);
    new(obj) TGswSampleFFT(params, all_samples);
}

//destructor content
EXPORT void destroy_TGswSampleFFT(TGswSampleFFT *obj) {
    int32_t k = obj->k;
    int32_t l = obj->l;
    delete_TLweSampleFFT_array((k + 1) * l, obj->all_samples);
    obj->~TGswSampleFFT();
}


// For all the kpl TLWE samples composing the TGSW sample
// It computes the inverse FFT of the coefficients of the TLWE sample
EXPORT void tGswToFFTConvert(TGswSampleFFT *result, const TGswSample *source, const TGswParams *params) {
    const int32_t kpl = params->kpl;

    for (int32_t p = 0; p < kpl; p++)
        tLweToFFTConvert(result->all_samples + p, source->all_sample + p, params->tlwe_params);
}

// For all the kpl TLWE samples composing the TGSW sample
// It computes the FFT of the coefficients of the TLWEfft sample
EXPORT void tGswFromFFTConvert(TGswSample *result, const TGswSampleFFT *source, const TGswParams *params) {
    const int32_t kpl = params->kpl;

    for (int32_t p = 0; p < kpl; p++)
        tLweFromFFTConvert(result->all_sample + p, source->all_samples + p, params->tlwe_params);
}



// result = result + H
EXPORT void tGswFFTAddH(TGswSampleFFT *result, const TGswParams *params) {
    const int32_t k = params->tlwe_params->k;
    const int32_t l = params->l;

    for (int32_t j = 0; j < l; j++) {
        Torus32 hj = params->h[j];
        for (int32_t i = 0; i <= k; i++)
            LagrangeHalfCPolynomialAddTorusConstant(&result->sample[i][j].a[i], hj);
    }

}

// result = list of TLWE (0,0)
EXPORT void tGswFFTClear(TGswSampleFFT *result, const TGswParams *params) {
    const int32_t kpl = params->kpl;

    for (int32_t p = 0; p < kpl; p++)
        tLweFFTClear(result->all_samples + p);
}

// External product (*): accum = gsw (*) accum
EXPORT void tGswFFTExternMulToTLwe(TLweSample *accum, const TGswSampleFFT *gsw, const TGswParams *params) {
    const TLweParams *tlwe_params = params->tlwe_params;
    const int32_t k = tlwe_params->k;
    const int32_t l = params->l;
    const int32_t kpl = params->kpl;

    //TODO attention, improve these new/delete...
    IntPolynomial *deca = new_IntPolynomial_array(kpl); //decomposed accumulator // Emerson Note: I think this should be kpl * k. tGswTorus32PolynomialDecompH implementation jumps by k
    LagrangeHalfCPolynomial *decaFFT = new_LagrangeHalfCPolynomial_array(kpl); //fft version
    TLweSampleFFT *tmpa = new_TLweSampleFFT(tlwe_params);

    // IntPolynomial_Collapsed deca = IntPolynomial_Collapsed[kpl];
    // LagrangeHalfCPolynomial_Collapsed decaFFT = LagrangeHalfCPolynomial_Collapsed[kpl];

    for (int32_t i = 0; i <= k; i++)
        tGswTorus32PolynomialDecompH(deca + i * l, accum->a + i);
    for (int32_t p = 0; p < kpl; p++)
        IntPolynomial_ifft(decaFFT[p].coefsC, deca[p].coefs);

    tLweFFTClear(tmpa);
    for (int32_t p = 0; p < kpl; p++) {
        tLweFFTAddMulRTo(tmpa, decaFFT + p, gsw->all_samples + p, tlwe_params); // TODO: Convert to HLS
    }
    tLweFromFFTConvert(accum, tmpa, tlwe_params); // TODO: Convert to HLS

    delete_TLweSampleFFT(tmpa);
    delete_LagrangeHalfCPolynomial_array(kpl, decaFFT);
    delete_IntPolynomial_array(kpl, deca);





    // cl::Buffer result_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(LagrangeHalfCPolynomial_Collapsed));
    // cl::Buffer p_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(IntPolynomial_Collapsed));

	// cplx *result_map = (cplx *)fpga.q.enqueueMapBuffer(result_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(LagrangeHalfCPolynomial_Collapsed));
	// int32_t *p_map = (int32_t *)fpga.q.enqueueMapBuffer(p_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(IntPolynomial_Collapsed));

	// memcpy(p_map, p, sizeof(int32_t) * Value_N);

	// fpga.k_IntPolynomial_ifft.setArg(0, result_buf);
	// fpga.k_IntPolynomial_ifft.setArg(1, p_buf);

	// fpga.q.enqueueMigrateMemObjects({ result_buf, p_buf }, 0 /* 0 means from host*/);
	// fpga.q.enqueueTask(fpga.k_IntPolynomial_ifft);
	// fpga.q.enqueueMigrateMemObjects({ result_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	// fpga.q.finish();

	// memcpy(result, result_map, sizeof(cplx) * Value_Ns2);
}

// result = (X^ai -1)*bki
/*
//This function is not used, but may become handy in a future release
//
EXPORT void tGswFFTMulByXaiMinusOne(TGswSampleFFT* result, const int32_t ai, const TGswSampleFFT* bki, const TGswParams* params) {
    const TLweParams* tlwe_params=params->tlwe_params;
    const int32_t k = tlwe_params->k;
    //const int32_t l = params->l;
    const int32_t kpl = params->kpl;
    const int32_t N = tlwe_params->N;
    //on calcule x^ai-1 en fft
    //TODO attention, this prevents parallelization...
    //TODO: parallelization
    static LagrangeHalfCPolynomial* xaim1=new_LagrangeHalfCPolynomial(N);
    LagrangeHalfCPolynomialSetXaiMinusOne(xaim1,ai);
    for (int32_t p=0; p<kpl; p++) {
        const LagrangeHalfCPolynomial* in_s = bki->all_samples[p].a;
        LagrangeHalfCPolynomial* out_s = result->all_samples[p].a;
        for (int32_t j=0; j<=k; j++)
            LagrangeHalfCPolynomialMul(&out_s[j], xaim1, &in_s[j]);
    }
}
*/

//-------------------------------------------------------------------------------------
// autogenerated memory-related functions
//-------------------------------------------------------------------------------------

USE_DEFAULT_CONSTRUCTOR_DESTRUCTOR_IMPLEMENTATIONS1(TGswSampleFFT, TGswParams);

//
//----------------------------------------------------------------------------------------







#if 0
// BOOTSTRAPPING (as in CGGI16b - algo 3)
//  - modswitch: torus coefs multiplied by N/2
//  - set the test vector
//  - blind rotation by the phase
//  - sample extract
//  - keyswitch
EXPORT void tfhe_bootstrapFFT(LweSample* result, const LweBootstrappingKeyFFT* bk, Torus32 mu1, Torus32 mu0, const LweSample* x){
    const Torus32 ab=(mu1+mu0)/2;
    const Torus32 aa = mu0-ab; // aa=(mu1-mu0)/2;
    const TGswParams* bk_params = bk->bk_params;
    const TLweParams* accum_params = bk_params->tlwe_params;
    const LweParams* extract_params = &accum_params->extracted_lweparams;
    const LweParams* in_out_params = bk->in_out_params;
    const int32_t n=in_out_params->n;
    const int32_t N=accum_params->N;
    const int32_t Ns2=N/2;
    const int32_t Nx2= 2*N;


    // Set the test vector (aa + aaX + ... + aaX^{N/2-1} -aaX^{N/2} - ... -aaX^{N-1})*X^{b}
    TorusPolynomial* testvect=new_TorusPolynomial();
    TorusPolynomial* testvectbis=new_TorusPolynomial();

    int32_t barb=modSwitchFromTorus32(x->b,Nx2);
    //je definis le test vector (multiplié par a inclus !
    for (int32_t i=0;i<Ns2;i++)
       testvect->coefsT[i]=aa;
    for (int32_t i=Ns2;i<N;i++)
       testvect->coefsT[i]=-aa;
    torusPolynomialMulByXai(testvectbis, barb, testvect);



    // Accumulateur acc = fft((0, testvect))
    TLweSample* acc = new_TLweSample(accum_params);

    // acc will be used for tfhe_bootstrapFFT, acc1=acc will be used for tfhe_bootstrap
    tLweNoiselessTrivial(acc, testvectbis, accum_params);

    TGswSample* temp = new_TGswSample(bk_params);
    TGswSampleFFT* tempFFT = new_TGswSampleFFT(bk_params);


    // Blind rotation
//NICOLAS: j'ai ajouté ce bloc
#ifndef NDEBUG
    TorusPolynomial* phase = new_TorusPolynomial();
    int32_t correctOffset = barb;
    cout << "starting the test..." << endl;
#endif
    // the index 1 is given when we don't use the fft
    for (int32_t i=0; i<n; i++) {
        int32_t bara=modSwitchFromTorus32(-x->a[i],Nx2);

        if (bara!=0) {
            tGswFFTMulByXaiMinusOne(tempFFT, bara, bk->bkFFT+i, bk_params);
            tGswFFTAddH(tempFFT, bk_params);
            tGswFFTExternMulToTLwe(acc, tempFFT, bk_params);
        }

//NICOLAS: et surtout, j'ai ajouté celui-ci!
#ifndef NDEBUG
    tLwePhase(phase,acc,debug_accum_key);  //celui-ci, c'est la phase de acc (FFT)
    if (debug_in_key->key[i]==1) correctOffset = (correctOffset+bara)%Nx2;
        torusPolynomialMulByXai(testvectbis, correctOffset, testvect); //celui-ci, c'est la phase idéale (calculée sans bruit avec la clé privée)
    for (int32_t j=0; j<N; j++) {
           printf("Iteration %d, index %d: phase %d vs noiseless %d\n",i,j,phase->coefsT[j], testvectbis->coefsT[j]);
    }
#endif

    }


    // Sample extract
    LweSample* u = new_LweSample(extract_params);
    tLweExtractLweSample(u, acc, extract_params, accum_params);
    u->b += ab;


    // KeySwitching
    lweKeySwitch(result, bk->ks, u);



    delete_LweSample(u);
    delete_TGswSampleFFT(tempFFT);
    delete_TGswSample(temp);
    delete_TLweSample(acc);
    delete_TorusPolynomial(testvectbis);
    delete_TorusPolynomial(testvect);
}
#endif

#undef INCLUDE_ALL
