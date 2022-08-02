#ifndef TFHE_CORE_H
#define TFHE_CORE_H

///@file
///@brief This file declares only the structures names
#include <stdint.h>

#define EXPORT extern "C"
#include "tfhe_generic_templates.h"
#include <complex>
typedef std::complex< double > cplx; // https://stackoverflow.com/a/31800404

constexpr static int64_t Value_N = 1024;
constexpr static int64_t Value_2N = 2 * Value_N;
constexpr static int64_t Value_Ns2 = Value_N / 2;
constexpr static int32_t Value_k = 1;
constexpr static int32_t Value_kpl = 3;
constexpr static int32_t Value_l = 3;
constexpr static int32_t Value_Bgbit = 7;
constexpr static int32_t Value_maskMod = 127;
constexpr static int32_t Value_halfBg = 64;
constexpr static uint32_t Value_offset = 2164391936; // tgsw.cpp - TGswParams::TGswParams()

EXPORT void die_dramatically(const char* message);


// Idea:
// we may want to represent an element x of the real torus by
// the integer rint(2^32.x) modulo 2^32
//  -- addition, subtraction and integer combinations are native operation
//  -- modulo 1 is mapped to mod 2^32, which is also native!
// This looks much better than using float/doubles, where modulo 1 is not
// natural at all.
typedef int32_t Torus32; //avant uint32_t
//typedef int64_t Torus64; //avant uint64_t

// "Collapsed" representation of datatypes that avoid using structs
// Must be defined identixcally to `fpga_constants.h`
typedef int32_t IntPolynomial_Collapsed[Value_N];
typedef Torus32 TorusPolynomial_Collapsed[Value_N];
typedef cplx LagrangeHalfCPolynomial_Collapsed[Value_Ns2];

typedef struct {
    LagrangeHalfCPolynomial_Collapsed a[Value_k + 1]; ///< array of length k+1: mask + right term
    // TODO: Reimplement `b` once needed...
    // LagrangeHalfCPolynomial_Collapsed b; ///< alias of a[k] to get the right term
    double current_variance; ///< avg variance of the sample
} TLweSampleFFT_FPGA;

typedef struct {
    TorusPolynomial_Collapsed a[Value_k + 1]; ///< array of length k+1: mask + right term
    // TODO: Reimplement `b` once needed...
    // TorusPolynomial *b; ///< alias of a[k] to get the right term
    double current_variance; ///< avg variance of the sample
} TLweSample_FPGA;

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
