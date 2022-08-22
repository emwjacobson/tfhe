#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

///@file
///@brief This file contains the declaration of polynomials structures

#include "tfhe_core.h"

/** This structure represents an integer polynomial modulo X^N+1 */
typedef struct IntPolynomial {
   int32_t coefs[Value_N];
} IntPolynomial;


/** This structure represents an torus polynomial modulo X^N+1 */
typedef struct TorusPolynomial {
   Torus32 coefsT[Value_N];
} TorusPolynomial;

/**
 * This structure is used for FFT operations, and is a representation
 * over C of a polynomial in R[X]/X^N+1
 * This type is meant to be specialized, and all implementations of the structure must be compatible
 * (reinterpret_cast) with this one. Namely, they should contain at most 2 pointers
 */
/**
 * structure that represents a real polynomial P mod X^N+1
 * as the N/2 complex numbers:
 * P(w), P(w^3), ..., P(w^(N-1))
 * where w is exp(i.pi/N)
 */
typedef struct LagrangeHalfCPolynomial
{
   cplx coefsC[Value_Ns2];
   // void* precomp; // This variable is likely not used, but kept for compatability
} LagrangeHalfCPolynomial;

typedef struct TLweSampleFFT_FPGA {
    LagrangeHalfCPolynomial a[Value_k + 1]; ///< array of length k+1: mask + right term
    // TODO: Reimplement `b` once needed...
    // LagrangeHalfCPolynomial b; ///< alias of a[k] to get the right term
    double current_variance; ///< avg variance of the sample
} TLweSampleFFT_FPGA;

typedef struct TLweSample_FPGA {
    TorusPolynomial a[Value_k + 1]; ///< array of length k+1: mask + right term
    // TODO: Reimplement `b` once needed...
    // TorusPolynomial *b; ///< alias of a[k] to get the right term
    double current_variance; ///< avg variance of the sample
} TLweSample_FPGA;

typedef struct TGswSampleFFT_FPGA {
    TLweSampleFFT_FPGA all_samples[(Value_k+1) * Value_l]; ///< TLweSample* all_sample; (k+1)l TLwe Sample
    // TODO: Reimplement `sample` when needed
    // TLweSampleFFT_FPGA *sample; ///< accès optionnel aux différents blocs de taille l. (optional access to the various blocks of size l.)
    //double current_variance;
    int32_t k;
    int32_t l;
} TGswSampleFFT_FPGA;

EXPORT void fft_transform_reverse(double *real, double *imag);
EXPORT void fft_transform(double *real, double *imag);

//allocate memory space for a IntPolynomial
EXPORT IntPolynomial* alloc_IntPolynomial();
EXPORT IntPolynomial* alloc_IntPolynomial_array(int32_t nbelts);

//free memory space for a IntPolynomial
EXPORT void free_IntPolynomial(IntPolynomial* ptr);
EXPORT void free_IntPolynomial_array(int32_t nbelts, IntPolynomial* ptr);

//initialize the IntPolynomial structure
//(equivalent of the C++ constructor)
EXPORT void init_IntPolynomial(IntPolynomial* obj);
EXPORT void init_IntPolynomial_array(int32_t nbelts, IntPolynomial* obj);

//destroys the IntPolynomial structure
//(equivalent of the C++ destructor)
EXPORT void destroy_IntPolynomial(IntPolynomial* obj);
EXPORT void destroy_IntPolynomial_array(int32_t nbelts, IntPolynomial* obj);

//allocates and initialize the IntPolynomial structure
//(equivalent of the C++ new)
EXPORT IntPolynomial* new_IntPolynomial();
EXPORT IntPolynomial* new_IntPolynomial_array(int32_t nbelts);

//destroys and frees the IntPolynomial structure
//(equivalent of the C++ delete)
EXPORT void delete_IntPolynomial(IntPolynomial* obj);
EXPORT void delete_IntPolynomial_array(int32_t nbelts, IntPolynomial* obj);

//allocate memory space for a TorusPolynomial
EXPORT TorusPolynomial* alloc_TorusPolynomial();
EXPORT TorusPolynomial* alloc_TorusPolynomial_array(int32_t nbelts);

//free memory space for a TorusPolynomial
EXPORT void free_TorusPolynomial(TorusPolynomial* ptr);
EXPORT void free_TorusPolynomial_array(int32_t nbelts, TorusPolynomial* ptr);

//initialize the TorusPolynomial structure
//(equivalent of the C++ constructor)
EXPORT void init_TorusPolynomial(TorusPolynomial* obj);
EXPORT void init_TorusPolynomial_array(int32_t nbelts, TorusPolynomial* obj);

//destroys the TorusPolynomial structure
//(equivalent of the C++ destructor)
EXPORT void destroy_TorusPolynomial(TorusPolynomial* obj);
EXPORT void destroy_TorusPolynomial_array(int32_t nbelts, TorusPolynomial* obj);

//allocates and initialize the TorusPolynomial structure
//(equivalent of the C++ new)
EXPORT TorusPolynomial* new_TorusPolynomial();
EXPORT TorusPolynomial* new_TorusPolynomial_array(int32_t nbelts);

//destroys and frees the TorusPolynomial structure
//(equivalent of the C++ delete)
EXPORT void delete_TorusPolynomial(TorusPolynomial* obj);
EXPORT void delete_TorusPolynomial_array(int32_t nbelts, TorusPolynomial* obj);

//allocate memory space for a LagrangeHalfCPolynomial
EXPORT LagrangeHalfCPolynomial* alloc_LagrangeHalfCPolynomial();
EXPORT LagrangeHalfCPolynomial* alloc_LagrangeHalfCPolynomial_array(int32_t nbelts);

//free memory space for a LagrangeHalfCPolynomial
EXPORT void free_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* ptr);
EXPORT void free_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* ptr);

//initialize the LagrangeHalfCPolynomial structure
//(equivalent of the C++ constructor)
EXPORT void init_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj);
EXPORT void init_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj);

//destroys the LagrangeHalfCPolynomial structure
//(equivalent of the C++ destructor)
EXPORT void destroy_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj);
EXPORT void destroy_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj);

//allocates and initialize the LagrangeHalfCPolynomial structure
//(equivalent of the C++ new)
EXPORT LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial();
EXPORT LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial_array(int32_t nbelts);

//destroys and frees the LagrangeHalfCPolynomial structure
//(equivalent of the C++ delete)
EXPORT void delete_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj);
EXPORT void delete_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj);

#endif //POLYNOMIALS_H
