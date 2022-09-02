#include <string.h>
#include "fpga_constants.h"

extern "C" {
  void tGswFFTExternMulToTLwe(TLweSample_FPGA *accum, const TGswSampleFFT_FPGA *gsw) {
    IntPolynomial deca[param_kpl];
    LagrangeHalfCPolynomial decaFFT[param_kpl];
    TLweSampleFFT_FPGA tmpa[param_kpl];
    TLweSampleFFT_FPGA tmpa_final;

    #pragma HLS array_partition variable=deca complete
    #pragma HLS array_partition variable=decaFFT complete
    #pragma HLS array_partition variable=tmpa complete

    loop_decomp: for(int i=0; i<=param_k; i++) {
      #pragma HLS unroll
      tGswTorus32PolynomialDecompH(&deca[i * param_l], &accum->a[i]);
    }
    loop_ifft: for(int p=0; p<param_kpl; p++) {
      #pragma HLS unroll
      IntPolynomial_ifft(&decaFFT[p], &deca[p]);
    }
    loop_clear: for(int p=0; p<param_kpl; p++) {
      #pragma HLS unroll
      tLweFFTClear(&tmpa[p]);
    }
    loop_rto: for(int p=0; p<param_kpl; p++) {
      #pragma HLS unroll
      tLweFFTAddMulRTo(&tmpa[p], &decaFFT[p], &gsw->all_samples[p]);
    }

    for(int p=0; p<param_kpl; p++) {
      for(int k=0; k<=param_k; k++) {
        #pragma HLS unroll
        for(int i=0; i<param_Ns2; i++) {
          #pragma HLS loop_flatten
          tmpa_final.a[k].coefsC[i] += tmpa[p].a[k].coefsC[i];
        }
      }
    }

    tLweFromFFTConvert(accum, &tmpa_final);
  }
}