#ifndef TFHE_CORE_H
#define TFHE_CORE_H

///@file
///@brief This file declares only the structures names
#include <stdint.h>

//not very important, but all the functions exported in the output library
//should use the C naming convention, for inter-compiler compatibility
//reasons, or to bind it with non C++ code.
#ifdef __cplusplus
#define EXPORT extern "C"
#include "tfhe_generic_templates.h"
#else
#define EXPORT
#endif

EXPORT void die_dramatically(const char* message);

// These are the parameter set provided in CGGI2019.
// Currently (in 2020), the security of these parameters is estimated to lambda = 129-bit security
// (w.r.t bkz-sieve model, + additional hybrid attack models)
#define PARAM_N 1024
#define PARAM_k 1
#define PARAM_n 630
#define PARAM_bk_l 3
#define PARAM_bk_Bgbit 7
#define PARAM_ks_basebit 2
#define PARAM_ks_basebit 2
#define PARAM_ks_length 8
#define PARAM_ks_stdev pow(2.,-15)
#define PARAM_bk_stdev pow(2.,-25)
#define PARAM_max_stdev 0.012467

// Idea:
// we may want to represent an element x of the real torus by
// the integer rint(2^32.x) modulo 2^32
//  -- addition, subtraction and integer combinations are native operation
//  -- modulo 1 is mapped to mod 2^32, which is also native!
// This looks much better than using float/doubles, where modulo 1 is not
// natural at all.
typedef int32_t Torus32; //avant uint32_t
//typedef int64_t Torus64; //avant uint64_t

// ------------------------------------------------------------------------------------
// Container structs
// ------------------------------------------------------------------------------------

typedef struct {
	Torus32 a[PARAM_n]; //-- the n coefs of the mask
	Torus32 b;  //
	double current_variance; //-- average noise of the sample
} LweSample_Container;

typedef struct {
	int32_t n;
	double alpha_min;
	double alpha_max;
} LweParams_Container;

typedef struct {
	int32_t N;
	int32_t k;
	double alpha_min;
	double alpha_max;
	LweParams_Container extracted_lweparams;
} TLweParams_Container;

typedef struct {
	int32_t l;
	int32_t Bgbit;
	int32_t Bg;
	int32_t halfBg;
	uint32_t maskMod;
	TLweParams_Container tlwe_params;
	int32_t kpl;
	Torus32 h[PARAM_bk_l];
	uint32_t offset;
} TGswParams_Container;

typedef struct {
	int32_t ks_t;
	int32_t ks_basebit;
	LweParams_Container in_out_params;
	TGswParams_Container tgsw_params;
} TFheGateBootstrappingParameterSet_Container;

// Was slightly confused making this one

typedef struct {
	int32_t n;
	int32_t t;
	int32_t basebit;
	int32_t base;
	LweParams_Container out_params;
	LweSample_Container ks0_raw[PARAM_N * PARAM_ks_length * (1 << PARAM_ks_basebit)];
	// TODO: Figure out what to do with the ks pointer arrays.
	LweSample* ks1_raw[PARAM_N]; // This is left as a pointer, as it points to data in ks0_raw
	LweSample** ks[PARAM_N]; // This is left as pointers to pointers as it points to ks1_raw, which in itself is an array of pointers
} LweKeySwitchKey_Container;

typedef struct {
	LweParams_Container in_out_params;
	TGswParams_Container bk_params;
	TLweParams_Container accum_params;
	LweParams_Container extract_params;
	TGswSample_Container bk[PARAM_n]; ///< the bootstrapping key (s->s")
	LweKeySwitchKey_Container ks; ///< the keyswitch key (s'->s)
} LweBootstrappingKey_Container;

typedef struct {
	TFheGateBootstrappingParameterSet_Container params;
  LweBootstrappingKey_Container bk;
  LweBootstrappingKeyFFT_Container bkFFT;
} TFheGateBootstrappingCloudKeySet_Container;

// ------------------------------------------------------------------------------------
// Container structs END
// ------------------------------------------------------------------------------------


struct LweParams;
struct LweKey;
struct LweSample;
struct LweKeySwitchKey;
struct TLweParams;
struct TLweKey;
struct TLweSample;
struct TLweSampleFFT;
struct TGswParams;
struct TGswKey;
struct TGswSample;
struct TGswSampleFFT;
struct LweBootstrappingKey;
struct LweBootstrappingKeyFFT;
struct IntPolynomial;
struct TorusPolynomial;
struct LagrangeHalfCPolynomial;
struct TFheGateBootstrappingParameterSet;
struct TFheGateBootstrappingCloudKeySet;
struct TFheGateBootstrappingSecretKeySet;

//this is for compatibility with C code, to be able to use
//"LweParams" as a type and not "struct LweParams"
typedef struct LweParams           LweParams;
typedef struct LweKey              LweKey;
typedef struct LweSample           LweSample;
typedef struct LweKeySwitchKey     LweKeySwitchKey;
typedef struct TLweParams       TLweParams;
typedef struct TLweKey          TLweKey;
typedef struct TLweSample       TLweSample;
typedef struct TLweSampleFFT       TLweSampleFFT;
typedef struct TGswParams       TGswParams;
typedef struct TGswKey          TGswKey;
typedef struct TGswSample       TGswSample;
typedef struct TGswSampleFFT       TGswSampleFFT;
typedef struct LweBootstrappingKey LweBootstrappingKey;
typedef struct LweBootstrappingKeyFFT LweBootstrappingKeyFFT;
typedef struct IntPolynomial	   IntPolynomial;
typedef struct TorusPolynomial	   TorusPolynomial;
typedef struct LagrangeHalfCPolynomial	   LagrangeHalfCPolynomial;
typedef struct TFheGateBootstrappingParameterSet TFheGateBootstrappingParameterSet;
typedef struct TFheGateBootstrappingCloudKeySet TFheGateBootstrappingCloudKeySet;
typedef struct TFheGateBootstrappingSecretKeySet TFheGateBootstrappingSecretKeySet;

#endif //TFHE_CORE_H
