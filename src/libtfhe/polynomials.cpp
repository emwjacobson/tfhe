#include <cassert>
#include <cmath>
#include <cstdlib>
#include "tfhe_core.h"
#include "polynomials_arithmetic.h"
#include "lagrangehalfc_arithmetic.h"
#include "numeric_functions.h"
#include "polynomials.h"
#include "fft.h"
#include "fpga.h"

using namespace std;

//allocate memory space for a LagrangeHalfCPolynomial
EXPORT LagrangeHalfCPolynomial* alloc_LagrangeHalfCPolynomial() {
    return (LagrangeHalfCPolynomial*) malloc(sizeof(LagrangeHalfCPolynomial));
}
EXPORT LagrangeHalfCPolynomial* alloc_LagrangeHalfCPolynomial_array(int32_t nbelts) {
    return (LagrangeHalfCPolynomial*) malloc(nbelts*sizeof(LagrangeHalfCPolynomial));
}

//free memory space for a LweKey
EXPORT void free_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* ptr) {
    free(ptr);
}
EXPORT void free_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* ptr) {
    free(ptr);
}

//initialize the key structure
//(equivalent of the C++ constructor)
EXPORT void init_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj) {
    new(obj) LagrangeHalfCPolynomial();
}
EXPORT void init_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj) {
    for (int32_t i=0; i<nbelts; i++) {
	new(obj+i) LagrangeHalfCPolynomial();
    }
}

//allocates and initialize the LagrangeHalfCPolynomial structure
//(equivalent of the C++ new)
EXPORT LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial() {
    LagrangeHalfCPolynomial* obj = alloc_LagrangeHalfCPolynomial();
    init_LagrangeHalfCPolynomial(obj);
    return obj;
}
EXPORT LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial_array(int32_t nbelts) {
    LagrangeHalfCPolynomial* obj = alloc_LagrangeHalfCPolynomial_array(nbelts);
    init_LagrangeHalfCPolynomial_array(nbelts,obj);
    return obj;
}

//destroys the LagrangeHalfCPolynomial structure
//(equivalent of the C++ destructor)
EXPORT void destroy_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj) {
    obj->~LagrangeHalfCPolynomial();
}
EXPORT void destroy_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj) {
    for (int32_t i=0; i<nbelts; i++) {
	(obj+i)->~LagrangeHalfCPolynomial();
    }
}

//destroys and frees the LagrangeHalfCPolynomial structure
//(equivalent of the C++ delete)
EXPORT void delete_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj) {
    destroy_LagrangeHalfCPolynomial(obj);
    free_LagrangeHalfCPolynomial(obj);
}
EXPORT void delete_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj) {
    destroy_LagrangeHalfCPolynomial_array(nbelts,obj);
    free_LagrangeHalfCPolynomial_array(nbelts,obj);
}


/** multiplication via direct FFT (it must know the implem of LagrangeHalfCPolynomial because of the tmp+1 notation */
EXPORT void torusPolynomialMultFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    const int32_t N = poly1->N;
    LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3);
    IntPolynomial_ifft(tmp+0,poly1);
    TorusPolynomial_ifft(tmp+1,poly2);
    LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
    TorusPolynomial_fft(result, tmp+2);
    delete_LagrangeHalfCPolynomial_array(3,tmp);
}
EXPORT void torusPolynomialAddMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    const int32_t N = poly1->N;
    LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3);
    TorusPolynomial* tmpr = new_TorusPolynomial(N);
    IntPolynomial_ifft(tmp+0,poly1);
    TorusPolynomial_ifft(tmp+1,poly2);
    LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
    TorusPolynomial_fft(tmpr, tmp+2);
    torusPolynomialAddTo(result, tmpr);
    delete_TorusPolynomial(tmpr);
    delete_LagrangeHalfCPolynomial_array(3,tmp);
}
EXPORT void torusPolynomialSubMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    const int32_t N = poly1->N;
    LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3);
    TorusPolynomial* tmpr = new_TorusPolynomial(N);
    IntPolynomial_ifft(tmp+0,poly1);
    TorusPolynomial_ifft(tmp+1,poly2);
    LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
    TorusPolynomial_fft(tmpr, tmp+2);
    torusPolynomialSubTo(result, tmpr);
    delete_TorusPolynomial(tmpr);
    delete_LagrangeHalfCPolynomial_array(3,tmp);
}

//MISC OPERATIONS
/** sets to zero */
EXPORT void LagrangeHalfCPolynomialClear(
	LagrangeHalfCPolynomial* reps) {
    for (int32_t i=0; i<Value_Ns2; i++)
	reps->coefsC[i] = 0;
}

EXPORT void LagrangeHalfCPolynomialSetTorusConstant(LagrangeHalfCPolynomial* result, const Torus32 mu) {
    cplx* b = result->coefsC;
    const cplx muc = t32tod(mu);
    for (int32_t j=0; j<Value_Ns2; j++)
    	b[j]=muc;
}

EXPORT void LagrangeHalfCPolynomialAddTorusConstant(LagrangeHalfCPolynomial* result, const Torus32 mu) {
    cplx* b = result->coefsC;
    const cplx muc = t32tod(mu);
    for (int32_t j=0; j<Value_Ns2; j++)
    	b[j]+=muc;
}

EXPORT void LagrangeHalfCPolynomialSetXaiMinusOne(LagrangeHalfCPolynomial* result, const int32_t ai) {
    const cplx* omegaxminus1 = fpga.omegaxminus1;
    for (int32_t i=0; i<Value_Ns2; i++)
	result->coefsC[i]=omegaxminus1[((2*i+1)*ai)%Value_2N];
}

/** termwise multiplication in Lagrange space */
EXPORT void LagrangeHalfCPolynomialMul(
	LagrangeHalfCPolynomial* result,
	const LagrangeHalfCPolynomial* a,
	const LagrangeHalfCPolynomial* b) {
    const cplx* aa = a->coefsC;
    const cplx* bb = b->coefsC;
    cplx* rr = result->coefsC;
    for (int32_t i=0; i<Value_Ns2; i++)
	rr[i] = aa[i]*bb[i];
}

/** termwise multiplication and addTo in Lagrange space */
EXPORT void LagrangeHalfCPolynomialAddMul(
	LagrangeHalfCPolynomial* accum,
	const LagrangeHalfCPolynomial* a,
	const LagrangeHalfCPolynomial* b)
{
    const cplx* aa = a->coefsC;
    const cplx* bb = b->coefsC;
    cplx* rr = accum->coefsC;
    for (int32_t i=0; i<Value_Ns2; i++)
	rr[i] += aa[i]*bb[i];
}

/** termwise multiplication and addTo in Lagrange space */
EXPORT void LagrangeHalfCPolynomialSubMul(
	LagrangeHalfCPolynomial* accum,
	const LagrangeHalfCPolynomial* a,
	const LagrangeHalfCPolynomial* b)
{
    const cplx* aa = a->coefsC;
    const cplx* bb = b->coefsC;
    cplx* rr = accum->coefsC;
    for (int32_t i=0; i<Value_Ns2; i++)
	rr[i] -= aa[i]*bb[i];
}

EXPORT void LagrangeHalfCPolynomialAddTo(
	LagrangeHalfCPolynomial* accum,
	const LagrangeHalfCPolynomial* a) {
    const cplx* aa = a->coefsC;
    cplx* rr = accum->coefsC;
    for (int32_t i=0; i<Value_Ns2; i++)
	rr[i] += aa[i];
}

EXPORT void IntPolynomial_ifft(LagrangeHalfCPolynomial* result, const IntPolynomial* p) {
    cplx* res = (cplx*)result->coefsC;
    const int32_t* a = p->coefs;

    double real_inout[Value_2N];
    double imag_inout[Value_2N];

    double* res_dbl=(double*) res;
    for (int32_t i=0; i<Value_N; i++) real_inout[i]=a[i]/2.;
    for (int32_t i=0; i<Value_N; i++) real_inout[Value_N+i]=-real_inout[i];
    for (int32_t i=0; i<Value_2N; i++) imag_inout[i]=0;

    fft_transform_reverse(real_inout, imag_inout);

    for (int32_t i=0; i<Value_N; i+=2) {
        res_dbl[i]=real_inout[i+1];
        res_dbl[i+1]=imag_inout[i+1];
    }
}

EXPORT void TorusPolynomial_ifft(LagrangeHalfCPolynomial* result, const TorusPolynomial* p) {
    cplx* res = result->coefsC;
    Torus32* a = p->coefsT;

    double real_inout[Value_2N];
    double imag_inout[Value_2N];

    static const double _2pm33 = 1./double(INT64_C(1)<<33);
    int32_t* aa = (int32_t*) a;
    for (int32_t i=0; i<Value_N; i++) real_inout[i]=aa[i]*_2pm33;
    for (int32_t i=0; i<Value_N; i++) real_inout[Value_N+i]=-real_inout[i];
    for (int32_t i=0; i<Value_2N; i++) imag_inout[i]=0;

    fft_transform_reverse(real_inout, imag_inout);

    for (int32_t i=0; i<Value_Ns2; i++) res[i]=cplx(real_inout[2*i+1],imag_inout[2*i+1]);
}

EXPORT void TorusPolynomial_fft(TorusPolynomial* result, const LagrangeHalfCPolynomial* p) {
    Torus32* res = result->coefsT;
    const cplx* a = p->coefsC;

    double real_inout[Value_2N];
    double imag_inout[Value_2N];

    static const double _2p32 = double(INT64_C(1)<<32);
    static const double _1sN = double(1)/double(Value_N);
    //double* a_dbl=(double*) a;
    for (int32_t i=0; i<Value_N; i++) real_inout[2*i]=0;
    for (int32_t i=0; i<Value_N; i++) imag_inout[2*i]=0;
    for (int32_t i=0; i<Value_Ns2; i++) real_inout[2*i+1]=real(a[i]);
    for (int32_t i=0; i<Value_Ns2; i++) imag_inout[2*i+1]=imag(a[i]);
    for (int32_t i=0; i<Value_Ns2; i++) real_inout[Value_2N-1-2*i]=real(a[i]);
    for (int32_t i=0; i<Value_Ns2; i++) imag_inout[Value_2N-1-2*i]=-imag(a[i]);
    fft_transform(real_inout, imag_inout);
    for (int32_t i=0; i<Value_N; i++) res[i]=Torus32(int64_t(real_inout[i]*_1sN*_2p32));
    //pas besoin du fmod... Torus32(int64_t(fmod(rev_out[i]*_1sN,1.)*_2p32));
}
