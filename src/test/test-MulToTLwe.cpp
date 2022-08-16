#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "tfhe.h"
#include "polynomials.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"
#include "fpga.h"
#include <cassert>
#include <CL/cl2.hpp>

using namespace std;

#define EPSILON 1e-10

//this function will be used to generate a fake TLWE sample
void tLweSampleUniform(TLweSample *result) {
  for(int i=0; i<=Value_k; i++)
    for (int32_t j = 0; j < Value_N; ++j)
      result->a[i].coefsT[j] = uniformTorus32_distrib(generator);;
  result->current_variance = rand() / double(RAND_MAX);
}

void tGswFFTExternMulToTLwe_CPU(TLweSample *accum, const TGswSampleFFT *gsw, const TGswParams *params) {
  const TLweParams *tlwe_params = params->tlwe_params;

  //TODO attention, improve these new/delete...
  IntPolynomial *deca = new_IntPolynomial_array(Value_kpl); //decomposed accumulator // Emerson Note: I think this should be kpl * k. tGswTorus32PolynomialDecompH implementation jumps by k
  LagrangeHalfCPolynomial *decaFFT = new_LagrangeHalfCPolynomial_array(Value_kpl); //fft version
  TLweSampleFFT *tmpa = new_TLweSampleFFT(tlwe_params);

  // IntPolynomial_Collapsed deca = IntPolynomial_Collapsed[kpl];
  // LagrangeHalfCPolynomial_Collapsed decaFFT = LagrangeHalfCPolynomial_Collapsed[kpl];

  for (int32_t i = 0; i <= Value_k; i++) {
      tGswTorus32PolynomialDecompH(deca + i * Value_l, accum->a + i);
  }
  for (int32_t p = 0; p < Value_kpl; p++)
      IntPolynomial_ifft(&decaFFT[p], &deca[p]);

  tLweFFTClear(tmpa);
  for (int32_t p = 0; p < Value_kpl; p++) {
      tLweFFTAddMulRTo(tmpa, decaFFT + p, gsw->all_samples + p); // TODO: Convert to HLS
  }
  tLweFromFFTConvert(accum, tmpa);

  delete_TLweSampleFFT(tmpa);
  delete_LagrangeHalfCPolynomial_array(Value_kpl, decaFFT);
  delete_IntPolynomial_array(Value_kpl, deca);
}
void tGswFFTExternMulToTLwe_FPGA(TLweSample *accum, const TGswSampleFFT *gsw) {
  // Create buffers and map to host memory
  cl::Buffer accum_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TLweSample_FPGA));
  cl::Buffer gsw_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TGswSampleFFT_FPGA));

	TLweSample_FPGA *accum_map = (TLweSample_FPGA *)fpga.q.enqueueMapBuffer(accum_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TLweSample_FPGA));
	TGswSampleFFT_FPGA *gsw_map = (TGswSampleFFT_FPGA *)fpga.q.enqueueMapBuffer(gsw_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TGswSampleFFT_FPGA));

  // Copy data from arguments to buffers
  accum_map->current_variance = accum->current_variance;
  for(int i=0; i<=Value_k; i++) {
      for(int j=0; j<Value_N; j++) {
          accum_map->a[i][j] = accum->a[i].coefsT[j];
      }
  }

  gsw_map->k = gsw->k;
  gsw_map->l = gsw->l;
  for(int i=0; i<(Value_k+1)*Value_l; i++) {
      gsw_map->all_samples[i].current_variance = gsw->all_samples[i].current_variance;
      for(int j=0; j<=Value_k; j++) {
          for(int k=0; k<Value_Ns2; k++) {
              gsw_map->all_samples[i].a[j][k] = gsw->all_samples[i].a[j].coefsC[k];
          }
      }
  }

    // Set parameters and launch kernel
	fpga.k_tGswFFTExternMulToTLwe.setArg(0, accum_buf);
	fpga.k_tGswFFTExternMulToTLwe.setArg(1, gsw_buf);

	fpga.q.enqueueMigrateMemObjects({ accum_buf, gsw_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_tGswFFTExternMulToTLwe);
	fpga.q.enqueueMigrateMemObjects({ accum_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

  // Copy data back
  accum->current_variance = accum_map->current_variance;
  for(int i=0; i<=Value_k; i++) {
      for(int j=0; j<Value_N; j++) {
          accum->a[i].coefsT[j] = accum_map->a[i][j];
      }
  }
}
void test_tGswFFTExternMulToTLwe() {
  printf("Starting tGswFFTExternMulToTLwe\n");

  int32_t minimum_lambda = 100;
  TFheGateBootstrappingParameterSet *params = new_default_gate_bootstrapping_parameters(minimum_lambda);
  // const LweParams *in_out_params = params->in_out_params;
  TFheGateBootstrappingSecretKeySet *keyset = new_random_gate_bootstrapping_secret_keyset(params);
  TFheGateBootstrappingCloudKeySet *cloud = (TFheGateBootstrappingCloudKeySet *)(&keyset->cloud);

  TLweSample *result_CPU = new_TLweSample(params->tgsw_params->tlwe_params);
  TLweSample *result_FPGA = new_TLweSample(params->tgsw_params->tlwe_params);

  tLweSampleUniform(result_CPU);
  result_FPGA->current_variance = result_CPU->current_variance;
  for(int i=0; i<=Value_k; i++) {
    torusPolynomialCopy(&result_FPGA->a[i], &result_CPU->a[i]);
  }

  for(int i=0; i<=Value_k; i++) {
    for(int j=0; j<Value_N; j++) {
      if (result_CPU->a[i].coefsT[j] != result_FPGA->a[i].coefsT[j])
        printf("%i, %i\n", result_CPU->a[i].coefsT[j], result_FPGA->a[i].coefsT[j]);
    }
  }

  tGswFFTExternMulToTLwe_CPU(result_CPU, cloud->bkFFT->bkFFT, params->tgsw_params);
  tGswFFTExternMulToTLwe_FPGA(result_FPGA, cloud->bkFFT->bkFFT);

  // Used to verify that numbers are initialized and are not just 0s
  // for(int i=0; i<=Value_k; i++) {
  //   for(int j=0; j<=Value_k; j++) {
  //     for(int k=0; k<Value_Ns2; k++)
  //       printf("%f\n", (cloud->bkFFT->bkFFT->all_samples[i].a[j].coefsC[k].real()));
  //   }
  // }

  if (result_CPU->current_variance != result_FPGA->current_variance) {
    printf("ERROR variance does not match! CPU %f, FPGA %f\n", result_CPU->current_variance, result_FPGA->current_variance);
    exit(1);
  }
  for(int i=0; i<=Value_k; i++) {
    for(int j=0; j<Value_Ns2; j++) {
      if (result_CPU->a[i].coefsT[j] != result_FPGA->a[i].coefsT[j]) {
        printf("ERROR calculations don't match `%i` CPU %i, FPGA %i\n", i, result_CPU->a[i].coefsT[j], result_FPGA->a[i].coefsT[j]);
        exit(1);
      }
    }
  }

  printf("End tGswFFTExternMulToTLwe\n");
}


void tGswTorus32PolynomialDecompH_FPGA(IntPolynomial *result, const TorusPolynomial *sample) {
  uint32_t *buf = (uint32_t *)sample;

  cl::Buffer result_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(IntPolynomial_Collapsed) * Value_l);
  cl::Buffer sample_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TorusPolynomial_Collapsed));

	IntPolynomial_Collapsed *result_map = (IntPolynomial_Collapsed *)fpga.q.enqueueMapBuffer(result_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(IntPolynomial_Collapsed) * Value_l);
	int32_t *sample_map = (int32_t *)fpga.q.enqueueMapBuffer(sample_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TorusPolynomial_Collapsed));

	memcpy(sample_map, sample, sizeof(int32_t) * Value_N);

	fpga.k_tGswTorus32PolynomialDecompH.setArg(0, result_buf);
	fpga.k_tGswTorus32PolynomialDecompH.setArg(1, sample_buf);

	fpga.q.enqueueMigrateMemObjects({ result_buf, sample_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_tGswTorus32PolynomialDecompH);
	fpga.q.enqueueMigrateMemObjects({ result_buf, sample_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

	memcpy(buf, sample_map, sizeof(int32_t) * Value_N);
	memcpy(result, result_map, sizeof(IntPolynomial_Collapsed) * Value_l);
}
void IntPolynomial_ifft_FPGA(LagrangeHalfCPolynomial* result, const IntPolynomial* p) {
  cl::Buffer result_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(LagrangeHalfCPolynomial_Collapsed));
  cl::Buffer p_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(IntPolynomial_Collapsed));

	cplx *result_map = (cplx *)fpga.q.enqueueMapBuffer(result_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(LagrangeHalfCPolynomial_Collapsed));
	int32_t *p_map = (int32_t *)fpga.q.enqueueMapBuffer(p_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(IntPolynomial_Collapsed));

	memcpy(p_map, p, sizeof(int32_t) * Value_N);

	fpga.k_IntPolynomial_ifft.setArg(0, result_buf);
	fpga.k_IntPolynomial_ifft.setArg(1, p_buf);

	fpga.q.enqueueMigrateMemObjects({ result_buf, p_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_IntPolynomial_ifft);
	fpga.q.enqueueMigrateMemObjects({ result_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

	memcpy(result->coefsC, result_map, sizeof(cplx) * Value_Ns2);
}
void tLweFFTClear_FPGA(TLweSampleFFT *result) {
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
void tLweFFTAddMulRTo_FPGA(TLweSampleFFT *result, const LagrangeHalfCPolynomial *p, const TLweSampleFFT *sample) {
  cl::Buffer result_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TLweSampleFFT_FPGA));
  cl::Buffer p_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(LagrangeHalfCPolynomial_Collapsed));
  cl::Buffer sample_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TLweSampleFFT_FPGA));

	TLweSampleFFT_FPGA *result_map = (TLweSampleFFT_FPGA *)fpga.q.enqueueMapBuffer(result_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TLweSampleFFT_FPGA));
  cplx *p_map = (cplx *)fpga.q.enqueueMapBuffer(p_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(LagrangeHalfCPolynomial_Collapsed));
  TLweSampleFFT_FPGA *sample_map = (TLweSampleFFT_FPGA *)fpga.q.enqueueMapBuffer(sample_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TLweSampleFFT_FPGA));

  result_map->current_variance = result->current_variance;
  for(int i=0; i<=Value_k; i++) {
      for(int j=0; j<Value_Ns2; j++) {
          result_map->a[i][j] = result->a[i].coefsC[j];
      }
  }

  memcpy(p_map, p->coefsC, sizeof(cplx) * Value_Ns2);

  sample_map->current_variance = sample->current_variance;
  for(int i=0; i<=Value_k; i++) {
      for(int j=0; j<Value_Ns2; j++) {
          sample_map->a[i][j] = sample->a[i].coefsC[j];
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
          result->a[i].coefsC[j] = result_map->a[i][j];
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

  result->current_variance = result_map->current_variance;
  for(int i=0; i<=Value_k; i++) {
    for(int j=0; j<Value_N; j++) {
      result->a[i].coefsT[j] = result_map->a[i][j];
    }

    // for(int j=0; j<Value_Ns2; j++) {
    //     source->a[i].coefsC[j] = source_map->a[i][j];
    // }
  }
}

void manual_tGswFFTExternMulToTLwe() {
  printf("Starting manual_tGswFFTExternMulToTLwe\n");

  int32_t minimum_lambda = 100;
  TFheGateBootstrappingParameterSet *params = new_default_gate_bootstrapping_parameters(minimum_lambda);
  // const LweParams *in_out_params = params->in_out_params;
  TFheGateBootstrappingSecretKeySet *keyset = new_random_gate_bootstrapping_secret_keyset(params);
  TFheGateBootstrappingCloudKeySet *cloud = (TFheGateBootstrappingCloudKeySet *)(&keyset->cloud);

  TLweSample *result_CPU = new_TLweSample(params->tgsw_params->tlwe_params);
  TLweSample *result_FPGA = new_TLweSample(params->tgsw_params->tlwe_params);

  tLweSampleUniform(result_CPU);
  result_FPGA->current_variance = result_CPU->current_variance;
  for(int i=0; i<=Value_k; i++) {
    torusPolynomialCopy(&result_FPGA->a[i], &result_CPU->a[i]);
  }

  for(int i=0; i<=Value_k; i++) {
    for(int j=0; j<Value_N; j++) {
      if (result_CPU->a[i].coefsT[j] != result_FPGA->a[i].coefsT[j])
        printf("%i, %i\n", result_CPU->a[i].coefsT[j], result_FPGA->a[i].coefsT[j]);
    }
  }

  const TLweParams *tlwe_params = params->tgsw_params->tlwe_params;
  const int32_t k = tlwe_params->k;
  const int32_t l = params->tgsw_params->l;
  const int32_t kpl = params->tgsw_params->kpl;


  IntPolynomial *deca_CPU = new_IntPolynomial_array(kpl); //decomposed accumulator
  LagrangeHalfCPolynomial *decaFFT_CPU = new_LagrangeHalfCPolynomial_array(kpl); //fft version
  TLweSampleFFT *tmpa_CPU = new_TLweSampleFFT(tlwe_params);
  IntPolynomial *deca_FPGA = new_IntPolynomial_array(kpl); //decomposed accumulator
  LagrangeHalfCPolynomial *decaFFT_FPGA = new_LagrangeHalfCPolynomial_array(kpl); //fft version
  TLweSampleFFT *tmpa_FPGA = new_TLweSampleFFT(tlwe_params);


  for (int32_t i = 0; i <= k; i++) {
    tGswTorus32PolynomialDecompH(deca_CPU + i * l, result_CPU->a + i);
    tGswTorus32PolynomialDecompH_FPGA(deca_FPGA + i * l, result_FPGA->a + i);
  }
  for(int i=0; i<kpl; i++) {
    for(int j=0; j<Value_N; j++) {
      if (deca_CPU[i].coefs[j] != deca_FPGA[i].coefs[j]) {
        printf("ERROR deca doesnt match! Interation [%i][%i] CPU %i, FPGA %i\n", i, j, deca_CPU[i].coefs[j], deca_FPGA[i].coefs[j]);
        exit(1);
      }
    }
  }

  for (int32_t p = 0; p < kpl; p++) {
    IntPolynomial_ifft(&decaFFT_CPU[p], &deca_CPU[p]);
    IntPolynomial_ifft_FPGA(&decaFFT_FPGA[p], &deca_FPGA[p]);
  }
  for(int i=0; i<kpl; i++) {
    for(int j=0; j<Value_Ns2; j++) {
      if (fabs(decaFFT_CPU[i].coefsC[j].real() - decaFFT_FPGA[i].coefsC[j].real()) > EPSILON || fabs(decaFFT_CPU[i].coefsC[j].imag() - decaFFT_FPGA[i].coefsC[j].imag()) > EPSILON) {
        printf("ERROR decaFFT doesnt match! Interation [%i][%i] CPU (%.15f, %.15f), FPGA (%.15f, %.15f)\n", i, j, decaFFT_CPU[i].coefsC[j].real(), decaFFT_CPU[i].coefsC[j].imag(), decaFFT_FPGA[i].coefsC[j].real(), decaFFT_FPGA[i].coefsC[j].imag());
        exit(1);
      }
    }
  }

  tLweFFTClear(tmpa_CPU);
  tLweFFTClear_FPGA(tmpa_FPGA);
  for(int i=0; i<=Value_k; i++) {
    for(int j=0; j<Value_Ns2; j++) {
      if (tmpa_CPU->a[i].coefsC[j].real() != 0 || tmpa_CPU->a[i].coefsC[j].imag() != 0 || fabs(tmpa_CPU->a[i].coefsC[j].real() - tmpa_FPGA->a[i].coefsC[j].real()) > EPSILON || fabs(tmpa_CPU->a[i].coefsC[j].imag() - tmpa_FPGA->a[i].coefsC[j].imag()) > EPSILON) {
        printf("ERROR clear doesnt match! Interation [%i][%i] CPU (%.15f, %.15f), FPGA (%.15f, %.15f)\n", i, j, decaFFT_CPU[i].coefsC[j].real(), decaFFT_CPU[i].coefsC[j].imag(), decaFFT_FPGA[i].coefsC[j].real(), decaFFT_FPGA[i].coefsC[j].imag());
        exit(1);
      }
    }
  }

  for (int32_t p = 0; p < kpl; p++) {
      tLweFFTAddMulRTo(tmpa_CPU, decaFFT_CPU + p, cloud->bkFFT->bkFFT->all_samples + p);
      tLweFFTAddMulRTo_FPGA(tmpa_FPGA, decaFFT_FPGA + p, cloud->bkFFT->bkFFT->all_samples + p);
  }
  for(int i=0; i<=Value_k; i++) {
    for(int j=0; j<Value_Ns2; j++) {
      if (fabs(tmpa_CPU->a[i].coefsC[j].real() - tmpa_FPGA->a[i].coefsC[j].real()) > EPSILON || fabs(tmpa_CPU->a[i].coefsC[j].imag() - tmpa_FPGA->a[i].coefsC[j].imag()) > EPSILON) {
        printf("ERROR AddMulR doesnt match! Interation [%i][%i] CPU (%.15f, %.15f), FPGA (%.15f, %.15f)\n", i, j, decaFFT_CPU[i].coefsC[j].real(), decaFFT_CPU[i].coefsC[j].imag(), decaFFT_FPGA[i].coefsC[j].real(), decaFFT_FPGA[i].coefsC[j].imag());
        exit(1);
      }
    }
  }
  tLweFromFFTConvert(result_CPU, tmpa_CPU);
  tLweFromFFTConvert_FPGA(result_FPGA, tmpa_FPGA);
  for(int i=0; i<=Value_k; i++) {
    for(int j=0; j<Value_N; j++) {
      if (result_CPU->a[i].coefsT[j] != result_FPGA->a[i].coefsT[j]) {
        printf("ERROR FFTConvert does not match! Interation [%i][%i] CPU %i, FPGA %i\n", i, j, result_CPU->a[i].coefsT[j], result_FPGA->a[i].coefsT[j]);
        exit(1);
      }
    }
  }



  delete_TLweSampleFFT(tmpa_CPU);
  delete_LagrangeHalfCPolynomial_array(kpl, decaFFT_CPU);
  delete_IntPolynomial_array(kpl, deca_CPU);
  delete_TLweSampleFFT(tmpa_FPGA);
  delete_LagrangeHalfCPolynomial_array(kpl, decaFFT_FPGA);
  delete_IntPolynomial_array(kpl, deca_FPGA);





  if (result_CPU->current_variance != result_FPGA->current_variance) {
    printf("ERROR variance does not match! CPU %f, FPGA %f\n", result_CPU->current_variance, result_FPGA->current_variance);
    exit(1);
  }
  for(int i=0; i<=Value_k; i++) {
    for(int j=0; j<Value_Ns2; j++) {
      if (result_CPU->a[i].coefsT[j] != result_FPGA->a[i].coefsT[j]) {
        printf("ERROR calculations don't match `%i` CPU %i, FPGA %i\n", i, result_CPU->a[i].coefsT[j], result_FPGA->a[i].coefsT[j]);
        exit(1);
      }
    }
  }

  printf("End tGswFFTExternMulToTLwe\n");
}


int32_t main(int32_t argc, char **argv) {
  for(int i=0; i<100000; i++) {
    // test_tGswFFTExternMulToTLwe();
    manual_tGswFFTExternMulToTLwe();
  }
}
