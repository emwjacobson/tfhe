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

#include "fpga.h"

using namespace std;
#else
#undef EXPORT
#define EXPORT
#endif

#define EPSILON 0

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
        tLweFromFFTConvert(result->all_sample + p, source->all_samples + p);
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



void tGswTorus32PolynomialDecompH_FPGA(IntPolynomial *result, const TorusPolynomial *sample) {
  uint32_t *buf = (uint32_t *)sample;

  cl::Buffer result_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(IntPolynomial) * Value_l);
  cl::Buffer sample_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TorusPolynomial));

	IntPolynomial *result_map = (IntPolynomial *)fpga.q.enqueueMapBuffer(result_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(IntPolynomial) * Value_l);
	int32_t *sample_map = (int32_t *)fpga.q.enqueueMapBuffer(sample_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TorusPolynomial));

	memcpy(sample_map, sample, sizeof(int32_t) * Value_N);

	fpga.k_tGswTorus32PolynomialDecompH.setArg(0, result_buf);
	fpga.k_tGswTorus32PolynomialDecompH.setArg(1, sample_buf);

	fpga.q.enqueueMigrateMemObjects({ result_buf, sample_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_tGswTorus32PolynomialDecompH);
	fpga.q.enqueueMigrateMemObjects({ result_buf, sample_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

	memcpy(buf, sample_map, sizeof(int32_t) * Value_N);
	memcpy(result, result_map, sizeof(IntPolynomial) * Value_l);
}
void tGswTorus32PolynomialDecompH_FPGA_loop(IntPolynomial *deca, const TLweSample *accum) {
    // uint32_t *buf = (uint32_t *)sample;

    cl::Buffer deca_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(IntPolynomial) * Value_kpl);
    cl::Buffer accum_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TLweSample_FPGA));

    IntPolynomial *deca_map = (IntPolynomial *)fpga.q.enqueueMapBuffer(deca_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(IntPolynomial) * Value_kpl);
    TLweSample_FPGA *accum_map = (TLweSample_FPGA *)fpga.q.enqueueMapBuffer(accum_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TLweSample_FPGA));

    for(int i=0; i<Value_kpl; i++) {
        for(int j=0; j<Value_N; j++) {
            deca_map[i].coefs[j] = deca[i].coefs[j];
        }
    }

    accum_map->current_variance = accum->current_variance;
    for(int i=0; i<=Value_k; i++) {
        for(int j=0; j<Value_N; j++) {
            accum_map->a[i].coefsT[j] = accum->a[i].coefsT[j];
        }
    }

    fpga.k_tGswTorus32PolynomialDecompH_loop.setArg(0, deca_buf);
    fpga.k_tGswTorus32PolynomialDecompH_loop.setArg(1, accum_buf);

    fpga.q.enqueueMigrateMemObjects({ deca_buf, accum_buf }, 0 /* 0 means from host*/);
    fpga.q.enqueueTask(fpga.k_tGswTorus32PolynomialDecompH_loop);
    fpga.q.enqueueMigrateMemObjects({ deca_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

    fpga.q.finish();

    for(int i=0; i<Value_kpl; i++) {
        for(int j=0; j<Value_N; j++) {
            deca[i].coefs[j] = deca_map[i].coefs[j];
        }
    }
}
void IntPolynomial_ifft_FPGA(LagrangeHalfCPolynomial* result, const IntPolynomial* p) {
  cl::Buffer result_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(LagrangeHalfCPolynomial));
  cl::Buffer p_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(IntPolynomial));

	cplx *result_map = (cplx *)fpga.q.enqueueMapBuffer(result_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(LagrangeHalfCPolynomial));
	int32_t *p_map = (int32_t *)fpga.q.enqueueMapBuffer(p_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(IntPolynomial));

	memcpy(p_map, p, sizeof(int32_t) * Value_N);

	fpga.k_IntPolynomial_ifft.setArg(0, result_buf);
	fpga.k_IntPolynomial_ifft.setArg(1, p_buf);

	fpga.q.enqueueMigrateMemObjects({ result_buf, p_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_IntPolynomial_ifft);
	fpga.q.enqueueMigrateMemObjects({ result_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

	memcpy(result->coefsC, result_map, sizeof(cplx) * Value_Ns2);
}
void IntPolynomial_ifft_FPGA_loop(LagrangeHalfCPolynomial *decaFFT, IntPolynomial *deca) {
    cl::Buffer decaFFT_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(LagrangeHalfCPolynomial) * Value_kpl);
    cl::Buffer deca_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(IntPolynomial) * Value_kpl);

	LagrangeHalfCPolynomial *decaFFT_map = (LagrangeHalfCPolynomial *)fpga.q.enqueueMapBuffer(decaFFT_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(LagrangeHalfCPolynomial) * Value_kpl);
	IntPolynomial *deca_map = (IntPolynomial *)fpga.q.enqueueMapBuffer(deca_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(IntPolynomial) * Value_kpl);

	// memcpy(p_map, p, sizeof(int32_t) * Value_N);
    for(int i=0; i<Value_kpl; i++) {
        for(int j=0; j<Value_Ns2; j++) {
            decaFFT_map[i].coefsC[j] = decaFFT[i].coefsC[j];
        }

        for(int j=0; j<Value_N; j++) {
            deca_map[i].coefs[j] = deca[i].coefs[j];
        }
    }

	fpga.k_IntPolynomial_ifft_loop.setArg(0, decaFFT_buf);
	fpga.k_IntPolynomial_ifft_loop.setArg(1, deca_buf);

	fpga.q.enqueueMigrateMemObjects({ decaFFT_buf, deca_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_IntPolynomial_ifft_loop);
	fpga.q.enqueueMigrateMemObjects({ decaFFT_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

    for(int i=0; i<Value_kpl; i++) {
        for(int j=0; j<Value_Ns2; j++) {
            decaFFT[i].coefsC[j] = decaFFT_map[i].coefsC[j];
        }
    }
}
void tLweFFTClear_FPGA(TLweSampleFFT *result) {
  cl::Buffer result_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TLweSampleFFT_FPGA));

	TLweSampleFFT_FPGA *result_map = (TLweSampleFFT_FPGA *)fpga.q.enqueueMapBuffer(result_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TLweSampleFFT_FPGA));

  result_map->current_variance = result->current_variance;
  for(int i=0; i<Value_k; i++) {
      for(int j=0; j<Value_Ns2; j++) {
          result_map->a[i].coefsC[j] = result->a[i].coefsC[j];
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
          result->a[i].coefsC[j] = result_map->a[i].coefsC[j];
      }
  }
}
void tLweFFTAddMulRTo_FPGA(TLweSampleFFT *result, const LagrangeHalfCPolynomial *p, const TLweSampleFFT *sample) {
    cl::Buffer result_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TLweSampleFFT_FPGA));
    cl::Buffer p_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(LagrangeHalfCPolynomial));
    cl::Buffer sample_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TLweSampleFFT_FPGA));

    TLweSampleFFT_FPGA *result_map = (TLweSampleFFT_FPGA *)fpga.q.enqueueMapBuffer(result_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TLweSampleFFT_FPGA));
    cplx *p_map = (cplx *)fpga.q.enqueueMapBuffer(p_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(LagrangeHalfCPolynomial));
    TLweSampleFFT_FPGA *sample_map = (TLweSampleFFT_FPGA *)fpga.q.enqueueMapBuffer(sample_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TLweSampleFFT_FPGA));

    result_map->current_variance = result->current_variance;
    for(int i=0; i<=Value_k; i++) {
        for(int j=0; j<Value_Ns2; j++) {
            result_map->a[i].coefsC[j] = result->a[i].coefsC[j];
        }
    }

    memcpy(p_map, p->coefsC, sizeof(cplx) * Value_Ns2);

    sample_map->current_variance = sample->current_variance;
    for(int i=0; i<=Value_k; i++) {
        for(int j=0; j<Value_Ns2; j++) {
            sample_map->a[i].coefsC[j] = sample->a[i].coefsC[j];
        }
    }

	fpga.k_tLweFFTAddMulRTo.setArg(0, result_buf);
	fpga.k_tLweFFTAddMulRTo.setArg(1, p_buf);
	fpga.k_tLweFFTAddMulRTo.setArg(2, sample_buf);

	fpga.q.enqueueMigrateMemObjects({ result_buf, p_buf, sample_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_tLweFFTAddMulRTo);
	fpga.q.enqueueMigrateMemObjects({ result_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

    result->current_variance = result_map->current_variance;
    for(int i=0; i<=Value_k; i++) {
        for(int j=0; j<Value_Ns2; j++) {
            result->a[i].coefsC[j] = result_map->a[i].coefsC[j];
        }
    }
}
void tLweFFTAddMulRTo_FPGA_loop(TLweSampleFFT *tmpa, LagrangeHalfCPolynomial *decaFFT, const TGswSampleFFT *gsw) {
    cl::Buffer tmpa_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TLweSampleFFT_FPGA));
    cl::Buffer decaFFT_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(LagrangeHalfCPolynomial) * Value_kpl);
    cl::Buffer gsw_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TGswSampleFFT_FPGA));

    TLweSampleFFT_FPGA *tmpa_map = (TLweSampleFFT_FPGA *)fpga.q.enqueueMapBuffer(tmpa_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TLweSampleFFT_FPGA));
    LagrangeHalfCPolynomial *decaFFT_map = (LagrangeHalfCPolynomial *)fpga.q.enqueueMapBuffer(decaFFT_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(LagrangeHalfCPolynomial) * Value_kpl);
    TGswSampleFFT_FPGA *gsw_map = (TGswSampleFFT_FPGA *)fpga.q.enqueueMapBuffer(gsw_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TGswSampleFFT_FPGA));

    tmpa_map->current_variance = tmpa->current_variance;
    for(int i=0; i<=Value_k; i++) {
        for(int j=0; j<Value_Ns2; j++) {
            tmpa_map->a[i].coefsC[j] = tmpa->a[i].coefsC[j];
        }
    }

    for(int i=0; i<Value_kpl; i++) {
        for(int j=0; j<Value_Ns2; j++) {
            decaFFT_map[i].coefsC[j] = decaFFT[i].coefsC[j];
        }
    }

    gsw_map->k = gsw->k;
    gsw_map->l = gsw->l;
    for(int i=0; i<(Value_k+1)*Value_l; i++) {
        for (int j=0; j<=Value_k; j++) {
            for(int k=0; k<Value_Ns2; k++) {
                gsw_map->all_samples[i].a[j].coefsC[k] = gsw->all_samples[i].a[j].coefsC[k];
            }
        }
    }


	fpga.k_tLweFFTAddMulRTo_loop.setArg(0, tmpa_buf);
	fpga.k_tLweFFTAddMulRTo_loop.setArg(1, decaFFT_buf);
	fpga.k_tLweFFTAddMulRTo_loop.setArg(2, gsw_buf);

	fpga.q.enqueueMigrateMemObjects({ tmpa_buf, decaFFT_buf, gsw_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_tLweFFTAddMulRTo_loop);
	fpga.q.enqueueMigrateMemObjects({ tmpa_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

    tmpa->current_variance = tmpa_map->current_variance;
    for(int i=0; i<=Value_k; i++) {
        for(int j=0; j<Value_Ns2; j++) {
            tmpa->a[i].coefsC[j] = tmpa_map->a[i].coefsC[j];
        }
    }
}
void tLweFromFFTConvert_FPGA(TLweSample *result, const TLweSampleFFT *source) {
  cl::Buffer result_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TLweSample_FPGA));
  cl::Buffer source_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TLweSampleFFT_FPGA));

	TLweSample_FPGA *result_map = (TLweSample_FPGA *)fpga.q.enqueueMapBuffer(result_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TLweSample_FPGA));
	TLweSampleFFT_FPGA *source_map = (TLweSampleFFT_FPGA *)fpga.q.enqueueMapBuffer(source_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TLweSampleFFT_FPGA));

  result_map->current_variance = result->current_variance;
  source_map->current_variance = source->current_variance;

  for(int i=0; i<=Value_k; i++) {
    for(int j=0; j<Value_N; j++) {
      result_map->a[i].coefsT[j] = result->a[i].coefsT[j];
    }

    for(int j=0; j<Value_Ns2; j++) {
      source_map->a[i].coefsC[j] = source->a[i].coefsC[j];
    }
  }


	fpga.k_tLweFromFFTConvert.setArg(0, result_buf);
	fpga.k_tLweFromFFTConvert.setArg(1, source_buf);

	fpga.q.enqueueMigrateMemObjects({ result_buf, source_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_tLweFromFFTConvert);
	fpga.q.enqueueMigrateMemObjects({ result_buf, source_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

  result->current_variance = result_map->current_variance;
  for(int i=0; i<=Value_k; i++) {
    for(int j=0; j<Value_N; j++) {
      result->a[i].coefsT[j] = result_map->a[i].coefsT[j];
    }

    // for(int j=0; j<Value_Ns2; j++) {
    //     source->a[i].coefsC[j] = source_map->a[i][j];
    // }
  }
}

// External product (*): accum = gsw (*) accum
EXPORT void tGswFFTExternMulToTLwe(TLweSample *accum, const TGswSampleFFT *gsw, const TGswParams *params) {
    const TLweParams *tlwe_params = params->tlwe_params;
    const int32_t k = tlwe_params->k;
    const int32_t l = params->l;
    const int32_t kpl = params->kpl;

    TLweSample *fpga_accum = new_TLweSample(tlwe_params);
    TLweSample *cpu_accum = new_TLweSample(tlwe_params);


    ////////////////////////

    cl::Buffer accum_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TLweSample_FPGA));
    cl::Buffer gsw_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TGswSampleFFT_FPGA));

    TLweSample_FPGA *accum_map = (TLweSample_FPGA *)fpga.q.enqueueMapBuffer(accum_buf, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, sizeof(TLweSample_FPGA));
    TGswSampleFFT_FPGA *gsw_map = (TGswSampleFFT_FPGA *)fpga.q.enqueueMapBuffer(gsw_buf, CL_TRUE, CL_MAP_READ | CL_MAP_WRITE, 0, sizeof(TGswSampleFFT_FPGA));

    accum_map->current_variance = accum->current_variance;
    for (int i=0; i<=Value_k; i++) {
        for(int j=0; j<Value_N; j++) {
            accum_map->a[i].coefsT[j] = accum->a[i].coefsT[j];
        }
    }

    gsw_map->k = gsw->k;
    gsw_map->l = gsw->l;
    for(int i=0; i<(Value_k+1) * Value_l; i++) {
        gsw_map->all_samples[i].current_variance = gsw->all_samples[i].current_variance;
        for (int j=0; j<=Value_k; j++) {
            for(int k=0; k<Value_Ns2; k++) {
                gsw_map->all_samples[i].a[j].coefsC[k] = gsw->all_samples[i].a[j].coefsC[k];
            }
        }
    }

    fpga.k_tGswFFTExternMulToTLwe.setArg(0, accum_buf);
    fpga.k_tGswFFTExternMulToTLwe.setArg(1, gsw_buf);

    fpga.q.enqueueMigrateMemObjects({ accum_buf, gsw_buf }, 0);
    fpga.q.enqueueTask(fpga.k_tGswFFTExternMulToTLwe);
    fpga.q.enqueueMigrateMemObjects({ accum_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

    fpga.q.finish();

    fpga_accum->current_variance = accum_map->current_variance;
    for (int i=0; i<=Value_k; i++) {
        for(int j=0; j<Value_N; j++) {
            fpga_accum->a[i].coefsT[j] = accum_map->a[i].coefsT[j];
        }
    }







    ////////////////////////

    if (kpl != 6) {
        printf("ERROR: kpl != 6. kpl = %i\n", kpl);
    }
    if (k != 1) {
        printf("ERROR: k != 1. k = %i\n", k);
    }

    //TODO attention, improve these new/delete...
    IntPolynomial *deca_CPU = new_IntPolynomial_array(kpl); //decomposed accumulator
    // IntPolynomial *deca_FPGA = new_IntPolynomial_array(kpl);
    LagrangeHalfCPolynomial *decaFFT_CPU = new_LagrangeHalfCPolynomial_array(kpl); //fft version
    // LagrangeHalfCPolynomial *decaFFT_FPGA = new_LagrangeHalfCPolynomial_array(kpl);
    TLweSampleFFT *tmpa_CPU = new_TLweSampleFFT(tlwe_params);
    // TLweSampleFFT *tmpa_FPGA = new_TLweSampleFFT(tlwe_params);

    for (int32_t i = 0; i <= k; i++) {
        tGswTorus32PolynomialDecompH(deca_CPU + i * l, accum->a + i);
        // tGswTorus32PolynomialDecompH_FPGA(deca_FPGA + i * l, accum->a + i);
    }
    // tGswTorus32PolynomialDecompH_FPGA_loop(deca_FPGA, accum);
    // for(int i=0; i<Value_kpl; i++) {
    //     for(int j=0; j<Value_N; j++) {
    //         if (deca_CPU[i].coefs[j] != deca_FPGA[i].coefs[j]) {
    //             printf("ERROR deca doesnt match! Interation [%i][%i] CPU %i, FPGA %i\n", i, j, deca_CPU[i].coefs[j], deca_FPGA[i].coefs[j]);
    //             exit(1);
    //         }
    //     }
    // }
    for (int32_t p = 0; p < kpl; p++) {
        IntPolynomial_ifft(&decaFFT_CPU[p], &deca_CPU[p]);
        // IntPolynomial_ifft_FPGA(&decaFFT_FPGA[p], &deca_FPGA[p]);
    }
    // IntPolynomial_ifft_FPGA_loop(decaFFT_FPGA, deca_FPGA);
    // for(int i=0; i<kpl; i++) {
    //     for(int j=0; j<Value_Ns2; j++) {
    //         if (fabs(decaFFT_CPU[i].coefsC[j].real() - decaFFT_FPGA[i].coefsC[j].real()) > EPSILON || fabs(decaFFT_CPU[i].coefsC[j].imag() - decaFFT_FPGA[i].coefsC[j].imag()) > EPSILON) {
    //             printf("ERROR decaFFT doesnt match! Interation [%i][%i] CPU (%.15f, %.15f), FPGA (%.15f, %.15f)\n", i, j, decaFFT_CPU[i].coefsC[j].real(), decaFFT_CPU[i].coefsC[j].imag(), decaFFT_FPGA[i].coefsC[j].real(), decaFFT_FPGA[i].coefsC[j].imag());
    //             exit(1);
    //         }
    //     }
    // }

    tLweFFTClear(tmpa_CPU);
    // tLweFFTClear_FPGA(tmpa_FPGA);
    // for(int i=0; i<=Value_k; i++) {
    //     for(int j=0; j<Value_Ns2; j++) {
    //         if (tmpa_CPU->a[i].coefsC[j].real() != 0 || tmpa_CPU->a[i].coefsC[j].imag() != 0 || fabs(tmpa_CPU->a[i].coefsC[j].real() - tmpa_FPGA->a[i].coefsC[j].real()) > EPSILON || fabs(tmpa_CPU->a[i].coefsC[j].imag() - tmpa_FPGA->a[i].coefsC[j].imag()) > EPSILON) {
    //             printf("ERROR clear doesnt match! Interation [%i][%i] CPU (%.15f, %.15f), FPGA (%.15f, %.15f)\n", i, j, decaFFT_CPU[i].coefsC[j].real(), decaFFT_CPU[i].coefsC[j].imag(), decaFFT_FPGA[i].coefsC[j].real(), decaFFT_FPGA[i].coefsC[j].imag());
    //             exit(1);
    //         }
    //     }
    // }

    for (int32_t p = 0; p < kpl; p++) {
        tLweFFTAddMulRTo(tmpa_CPU, decaFFT_CPU + p, gsw->all_samples + p);
        // tLweFFTAddMulRTo_FPGA(tmpa_FPGA, decaFFT_FPGA + p, gsw->all_samples + p);
    }
    // tLweFFTAddMulRTo_FPGA_loop(tmpa_FPGA, decaFFT_FPGA, gsw);
    // for(int i=0; i<=Value_k; i++) {
    //     for(int j=0; j<Value_Ns2; j++) {
    //         if (fabs(tmpa_CPU->a[i].coefsC[j].real() - tmpa_FPGA->a[i].coefsC[j].real()) > EPSILON || fabs(tmpa_CPU->a[i].coefsC[j].imag() - tmpa_FPGA->a[i].coefsC[j].imag()) > EPSILON) {
    //             printf("ERROR AddMulR doesnt match! Interation [%i][%i] CPU (%.20f, %.20f), FPGA (%.20f, %.20f)\n", i, j, tmpa_CPU->a[i].coefsC[j].real(), tmpa_CPU->a[i].coefsC[j].imag(), tmpa_FPGA->a[i].coefsC[j].real(), tmpa_FPGA->a[i].coefsC[j].imag());
    //             exit(1);
    //         }
    //     }
    // }
    tLweFromFFTConvert(cpu_accum, tmpa_CPU);
    // tLweFromFFTConvert_FPGA(fpga_accum, tmpa_FPGA);
    // for(int i=0; i<=Value_k; i++) {
    //     for(int j=0; j<Value_N; j++) {
    //         if (cpu_accum->a[i].coefsT[j] != fpga_accum->a[i].coefsT[j]) {
    //             printf("ERROR FFTConvert does not match! Interation [%i][%i] CPU %i, FPGA %i\n", i, j, cpu_accum->a[i].coefsT[j], fpga_accum->a[i].coefsT[j]);
    //             exit(1);
    //         }
    //     }
    // }

    union d_convert {
        double dVal;
        uint64_t iVal;
    };

    d_convert c_cpu;
    c_cpu.dVal = cpu_accum->current_variance;
    d_convert c_fpga;
    c_fpga.dVal = fpga_accum->current_variance;

    // Compare
    for(int i=0; i<=Value_k; i++) {
        for(int j=0; j<Value_N; j++) {
            if (cpu_accum->current_variance != fpga_accum->current_variance) {
                printf("Error in variance! CPU: %.30f, FPGA: %.30f\n", cpu_accum->current_variance, fpga_accum->current_variance);
                printf("Error in variance! CPU: %lX, FPGA: %lX\n", c_cpu.iVal, c_fpga.iVal);
                exit(1);
            }
            if (cpu_accum->a[i].coefsT[j] != fpga_accum->a[i].coefsT[j]) {
                printf("BAD CALC [%i][%i] CPU: %i, FPGA: %i\n", i, j, cpu_accum->a[i].coefsT[j], fpga_accum->a[i].coefsT[j]);
                exit(1);
            }
        }
    }

    // Copy either cpu or fpga version (SHOULD be identical...)
    // tLweCopy(accum, cpu_accum, tlwe_params);
    tLweCopy(accum, fpga_accum, tlwe_params);

    delete_TLweSampleFFT(tmpa_CPU);
    // delete_TLweSampleFFT(tmpa_FPGA);
    delete_LagrangeHalfCPolynomial_array(kpl, decaFFT_CPU);
    // delete_LagrangeHalfCPolynomial_array(kpl, decaFFT_FPGA);
    delete_IntPolynomial_array(kpl, deca_CPU);
    // delete_IntPolynomial_array(kpl, deca_FPGA);
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
