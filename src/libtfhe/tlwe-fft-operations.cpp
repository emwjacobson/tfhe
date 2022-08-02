#ifndef TFHE_TEST_ENVIRONMENT
/* ***************************************************
   TLWE fft operations
 *************************************************** */

#include <random>
#include <cassert>
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

using namespace std;
#define INCLUDE_ALL

#else
#undef EXPORT
#define EXPORT
#endif


#include "fpga.h"


#if defined INCLUDE_ALL || defined INCLUDE_INIT_TLWESAMPLE_FFT
#undef INCLUDE_INIT_TLWESAMPLE_FFT
EXPORT void init_TLweSampleFFT(TLweSampleFFT *obj, const TLweParams *params) {
    //a is a table of k+1 polynomials, b is an alias for &a[k]
    const int32_t k = params->k;
    LagrangeHalfCPolynomial *a = new_LagrangeHalfCPolynomial_array(k + 1);
    double current_variance = 0;
    new(obj) TLweSampleFFT(params, a, current_variance);
}
#endif

#if defined INCLUDE_ALL || defined INCLUDE_DESTROY_TLWESAMPLE_FFT
#undef INCLUDE_DESTROY_TLWESAMPLE_FFT
EXPORT void destroy_TLweSampleFFT(TLweSampleFFT *obj) {
    const int32_t k = obj->k;
    delete_LagrangeHalfCPolynomial_array(k + 1, obj->a);
    obj->~TLweSampleFFT();
}
#endif


#if defined INCLUDE_ALL || defined INCLUDE_TLWE_TO_FFT_CONVERT
#undef INCLUDE_TLWE_TO_FFT_CONVERT
// Computes the inverse FFT of the coefficients of the TLWE sample
EXPORT void tLweToFFTConvert(TLweSampleFFT *result, const TLweSample *source, const TLweParams *params) {
    const int32_t k = params->k;

    for (int32_t i = 0; i <= k; ++i)
        TorusPolynomial_ifft(result->a[i].coefsC, source->a[i].coefsT);
    result->current_variance = source->current_variance;
}
#endif


#if defined INCLUDE_ALL || defined INCLUDE_TLWE_FROM_FFT_CONVERT
#undef INCLUDE_TLWE_FROM_FFT_CONVERT
// Computes the FFT of the coefficients of the TLWEfft sample
EXPORT void tLweFromFFTConvert(TLweSample *result, const TLweSampleFFT *source) {
    cl::Buffer result_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TLweSample_FPGA));
    cl::Buffer source_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TLweSampleFFT_FPGA));

	TLweSample_FPGA *result_map = (TLweSample_FPGA *)fpga.q.enqueueMapBuffer(result_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TLweSample_FPGA));
	TLweSampleFFT_FPGA *source_map = (TLweSampleFFT_FPGA *)fpga.q.enqueueMapBuffer(source_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TLweSampleFFT_FPGA));

    result_map->current_variance = result->current_variance;
    source_map->current_variance = source->current_variance;

    for(int i=0; i<=Value_k; i++) {
        for(int j=0; j<Value_N; j++) {
            result_map->a[i][j] = result->a[i].coefsT[j];
        }

        for(int j=0; j<Value_Ns2; j++) {
            source_map->a[i][j] = source->a[i].coefsC[j];
        }
    }


	fpga.k_tLweFromFFTConvert.setArg(0, result_buf);
	fpga.k_tLweFromFFTConvert.setArg(1, source_buf);

	fpga.q.enqueueMigrateMemObjects({ result_buf, source_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_tLweFromFFTConvert);
	fpga.q.enqueueMigrateMemObjects({ result_buf, source_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

    for(int i=0; i<=Value_k; i++) {
        for(int j=0; j<Value_N; j++) {
            result->a[i].coefsT[j] = result_map->a[i][j];
        }

        for(int j=0; j<Value_Ns2; j++) {
            source->a[i].coefsC[j] = source_map->a[i][j];
        }
    }
}
#endif


#if defined INCLUDE_ALL || defined INCLUDE_TLWE_FFT_CLEAR
#undef INCLUDE_TLWE_FFT_CLEAR
//Arithmetic operations on TLwe samples
/** result = (0,0) */
EXPORT void tLweFFTClear(TLweSampleFFT *result) {
    cl::Buffer result_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TLweSampleFFT_FPGA));

	TLweSampleFFT_FPGA *result_map = (TLweSampleFFT_FPGA *)fpga.q.enqueueMapBuffer(result_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TLweSampleFFT_FPGA));

    result_map->current_variance = result->current_variance;
    for(int i=0; i<Value_k; i++) {
        for(int j=0; j<Value_Ns2; j++) {
            result_map->a[i][j] = result->a[i].coefsC[j];
        }
    }

	fpga.k_tLweFFTClear.setArg(0, result_buf);

	fpga.q.enqueueMigrateMemObjects({ result_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_tLweFFTClear);
	fpga.q.enqueueMigrateMemObjects({ result_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

    result->current_variance = result_map->current_variance;
    for(int i=0; i<=Value_k; i++) {
        for(int j=0; j<Value_Ns2; j++) {
            result->a[i].coefsC[j] = result_map->a[i][j];
        }
    }
}
#endif


#if defined INCLUDE_ALL || defined INCLUDE_TLWE_FFT_ADDMULRTO
#undef INCLUDE_TLWE_FFT_ADDMULRTO
// result = result + p*sample
EXPORT void tLweFFTAddMulRTo(TLweSampleFFT *result, const LagrangeHalfCPolynomial *p, const TLweSampleFFT *sample,
                             const TLweParams *params) {
    const int32_t k = params->k;

    for (int32_t i = 0; i <= k; i++)
        LagrangeHalfCPolynomialAddMul(result->a + i, p, sample->a + i);
    //result->current_variance += sample->current_variance;
    //TODO: how to compute the variance correctly?
}
#endif


//autogenerated memory functions (they will always be included, even in
//tests)

USE_DEFAULT_CONSTRUCTOR_DESTRUCTOR_IMPLEMENTATIONS1(TLweSampleFFT, TLweParams);

#undef INCLUDE_ALL
