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

#define EPSILON 1e-5

//this function will be used to generate a fake TLWE sample
void tLweSampleUniform(TLweSample *result) {
  for(int i=0; i<=Value_k; i++)
    for (int32_t j = 0; j < Value_N; ++j)
      result->a[i].coefsT[j] = uniformTorus32_distrib(generator);;
  result->current_variance = rand() / double(RAND_MAX);
}

void tGswFFTExternMulToTLwe_CPU(TLweSample *accum, const TGswSampleFFT *gsw, const TGswParams *params) {
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

  for (int32_t i = 0; i <= k; i++) {
      tGswTorus32PolynomialDecompH(deca + i * l, accum->a + i);
  }
  for (int32_t p = 0; p < kpl; p++)
      IntPolynomial_ifft(&decaFFT[p], &deca[p]);

  tLweFFTClear(tmpa);
  for (int32_t p = 0; p < kpl; p++) {
      tLweFFTAddMulRTo(tmpa, decaFFT + p, gsw->all_samples + p); // TODO: Convert to HLS
  }
  tLweFromFFTConvert(accum, tmpa);

  delete_TLweSampleFFT(tmpa);
  delete_LagrangeHalfCPolynomial_array(kpl, decaFFT);
  delete_IntPolynomial_array(kpl, deca);
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

int32_t main(int32_t argc, char **argv) {
#ifndef NDEBUG
  cout << "DEBUG MODE!" << endl;
#endif

/*
  TODO:
    tGswFFTExternMulToTLwe - CPU vs FPGA
*/

  for(int i=0; i<100000; i++) {
    test_tGswFFTExternMulToTLwe();
  }

  return 0;
}
