#ifndef POLYNOMIALS_H
#define POLYNOMIALS_H

///@file
///@brief This file contains the declaration of polynomials structures

#include "tfhe_core.h"
#include "fft.h"

typedef struct {
   const int64_t N;
   const int64_t _2N;
   const int64_t Ns2;
} N_Values_t;

static const N_Values_t N_Values = { 1024, 2*1024, 1024/2 };

class FFT_Processor_nayuki {
    public:
    const int32_t _2N;
    const int32_t N;
    const int32_t Ns2;
    private:
    double* real_inout;
    double* imag_inout;
    void* tables_direct;
    void* tables_reverse;
    public:
    cplx* omegaxminus1;

    FFT_Processor_nayuki(const int32_t N);
    void execute_direct_torus32(Torus32* res, const cplx* a);
    ~FFT_Processor_nayuki();
};

/** This structure represents an integer polynomial modulo X^N+1 */
struct IntPolynomial {
   const int32_t N;
   int32_t* coefs;

#ifdef __cplusplus
   IntPolynomial(const int32_t N);
   ~IntPolynomial();
   IntPolynomial(const IntPolynomial&) = delete; //forbidden
   IntPolynomial* operator=(const IntPolynomial&) = delete; //forbidden
#endif
};


/** This structure represents an torus polynomial modulo X^N+1 */
struct TorusPolynomial {
   const int32_t N;
   Torus32* coefsT;

#ifdef __cplusplus
   TorusPolynomial(const int32_t N);
   ~TorusPolynomial();
   TorusPolynomial(const TorusPolynomial&) = delete; //forbidden
   TorusPolynomial* operator=(const TorusPolynomial&) = delete; //forbidden
#endif
};

extern FFT_Processor_nayuki fp1024_nayuki;

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
struct LagrangeHalfCPolynomial
{
   cplx* coefsC;
   FFT_Processor_nayuki* proc;

   LagrangeHalfCPolynomial(int32_t N);
   ~LagrangeHalfCPolynomial();
};

//allocate memory space for a IntPolynomial
EXPORT IntPolynomial* alloc_IntPolynomial();
EXPORT IntPolynomial* alloc_IntPolynomial_array(int32_t nbelts);

//free memory space for a IntPolynomial
EXPORT void free_IntPolynomial(IntPolynomial* ptr);
EXPORT void free_IntPolynomial_array(int32_t nbelts, IntPolynomial* ptr);

//initialize the IntPolynomial structure
//(equivalent of the C++ constructor)
EXPORT void init_IntPolynomial(IntPolynomial* obj, const int32_t N);
EXPORT void init_IntPolynomial_array(int32_t nbelts, IntPolynomial* obj, const int32_t N);

//destroys the IntPolynomial structure
//(equivalent of the C++ destructor)
EXPORT void destroy_IntPolynomial(IntPolynomial* obj);
EXPORT void destroy_IntPolynomial_array(int32_t nbelts, IntPolynomial* obj);

//allocates and initialize the IntPolynomial structure
//(equivalent of the C++ new)
EXPORT IntPolynomial* new_IntPolynomial(const int32_t N);
EXPORT IntPolynomial* new_IntPolynomial_array(int32_t nbelts, const int32_t N);

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
EXPORT void init_TorusPolynomial(TorusPolynomial* obj, const int32_t N);
EXPORT void init_TorusPolynomial_array(int32_t nbelts, TorusPolynomial* obj, const int32_t N);

//destroys the TorusPolynomial structure
//(equivalent of the C++ destructor)
EXPORT void destroy_TorusPolynomial(TorusPolynomial* obj);
EXPORT void destroy_TorusPolynomial_array(int32_t nbelts, TorusPolynomial* obj);

//allocates and initialize the TorusPolynomial structure
//(equivalent of the C++ new)
EXPORT TorusPolynomial* new_TorusPolynomial(const int32_t N);
EXPORT TorusPolynomial* new_TorusPolynomial_array(int32_t nbelts, const int32_t N);

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
EXPORT void init_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj, const int32_t N);
EXPORT void init_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj, const int32_t N);

//destroys the LagrangeHalfCPolynomial structure
//(equivalent of the C++ destructor)
EXPORT void destroy_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj);
EXPORT void destroy_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj);

//allocates and initialize the LagrangeHalfCPolynomial structure
//(equivalent of the C++ new)
EXPORT LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial(const int32_t N);
EXPORT LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial_array(int32_t nbelts, const int32_t N);

//destroys and frees the LagrangeHalfCPolynomial structure
//(equivalent of the C++ delete)
EXPORT void delete_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj);
EXPORT void delete_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj);

#endif //POLYNOMIALS_H
