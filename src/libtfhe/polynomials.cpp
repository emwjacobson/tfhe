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

LagrangeHalfCPolynomial::LagrangeHalfCPolynomial(const int32_t N) {
    assert(N==1024);
    coefsC = new cplx[N/2];
}

LagrangeHalfCPolynomial::~LagrangeHalfCPolynomial() {
    delete[] coefsC;
}

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
EXPORT void init_LagrangeHalfCPolynomial(LagrangeHalfCPolynomial* obj, const int32_t N) {
    new(obj) LagrangeHalfCPolynomial(N);
}
EXPORT void init_LagrangeHalfCPolynomial_array(int32_t nbelts, LagrangeHalfCPolynomial* obj, const int32_t N) {
    for (int32_t i=0; i<nbelts; i++) {
	new(obj+i) LagrangeHalfCPolynomial(N);
    }
}

//allocates and initialize the LagrangeHalfCPolynomial structure
//(equivalent of the C++ new)
EXPORT LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial(const int32_t N) {
    LagrangeHalfCPolynomial* obj = alloc_LagrangeHalfCPolynomial();
    init_LagrangeHalfCPolynomial(obj,N);
    return obj;
}
EXPORT LagrangeHalfCPolynomial* new_LagrangeHalfCPolynomial_array(int32_t nbelts, const int32_t N) {
    LagrangeHalfCPolynomial* obj = alloc_LagrangeHalfCPolynomial_array(nbelts);
    init_LagrangeHalfCPolynomial_array(nbelts,obj,N);
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
    LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3,N);
    IntPolynomial_ifft(tmp+0,poly1);
    TorusPolynomial_ifft(tmp+1,poly2);
    LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
    TorusPolynomial_fft(result, tmp+2);
    delete_LagrangeHalfCPolynomial_array(3,tmp);
}
EXPORT void torusPolynomialAddMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    const int32_t N = poly1->N;
    LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3,N);
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
    LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3,N);
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
    for (int32_t i=0; i<N_Values.Ns2; i++)
	reps->coefsC[i] = 0;
}

EXPORT void LagrangeHalfCPolynomialSetTorusConstant(LagrangeHalfCPolynomial* result, const Torus32 mu) {
    cplx* b = result->coefsC;
    const cplx muc = t32tod(mu);
    for (int32_t j=0; j<N_Values.Ns2; j++)
    	b[j]=muc;
}

EXPORT void LagrangeHalfCPolynomialAddTorusConstant(LagrangeHalfCPolynomial* result, const Torus32 mu) {
    cplx* b = result->coefsC;
    const cplx muc = t32tod(mu);
    for (int32_t j=0; j<N_Values.Ns2; j++)
    	b[j]+=muc;
}

EXPORT void LagrangeHalfCPolynomialSetXaiMinusOne(LagrangeHalfCPolynomial* result, const int32_t ai) {
    const cplx* omegaxminus1 = fpga.omegaxminus1;
    for (int32_t i=0; i<N_Values.Ns2; i++)
	result->coefsC[i]=omegaxminus1[((2*i+1)*ai)%N_Values._2N];
}

/** termwise multiplication in Lagrange space */
EXPORT void LagrangeHalfCPolynomialMul(
	LagrangeHalfCPolynomial* result,
	const LagrangeHalfCPolynomial* a,
	const LagrangeHalfCPolynomial* b) {
    cplx* aa = a->coefsC;
    cplx* bb = b->coefsC;
    cplx* rr = result->coefsC;
    for (int32_t i=0; i<N_Values.Ns2; i++)
	rr[i] = aa[i]*bb[i];
}

/** termwise multiplication and addTo in Lagrange space */
EXPORT void LagrangeHalfCPolynomialAddMul(
	LagrangeHalfCPolynomial* accum,
	const LagrangeHalfCPolynomial* a,
	const LagrangeHalfCPolynomial* b)
{
    cplx* aa = a->coefsC;
    cplx* bb = b->coefsC;
    cplx* rr = accum->coefsC;
    for (int32_t i=0; i<N_Values.Ns2; i++)
	rr[i] += aa[i]*bb[i];
}

/** termwise multiplication and addTo in Lagrange space */
EXPORT void LagrangeHalfCPolynomialSubMul(
	LagrangeHalfCPolynomial* accum,
	const LagrangeHalfCPolynomial* a,
	const LagrangeHalfCPolynomial* b)
{
    cplx* aa = a->coefsC;
    cplx* bb = b->coefsC;
    cplx* rr = accum->coefsC;
    for (int32_t i=0; i<N_Values.Ns2; i++)
	rr[i] -= aa[i]*bb[i];
}

EXPORT void LagrangeHalfCPolynomialAddTo(
	LagrangeHalfCPolynomial* accum,
	const LagrangeHalfCPolynomial* a) {
    cplx* aa = a->coefsC;
    cplx* rr = accum->coefsC;
    for (int32_t i=0; i<N_Values.Ns2; i++)
	rr[i] += aa[i];
}

void check_alternate_real(const double *real_inout, const double *imag_inout) {
#ifndef NDEBUG
    for (int32_t i=0; i<N_Values._2N; i++) assert(fabs(imag_inout[i])<1e-8);
    for (int32_t i=0; i<N_Values.N; i++) assert(fabs(real_inout[i]+real_inout[N_Values.N+i])<1e-9);
#endif
}

void check_conjugate_cplx(const double *real_inout, const double *imag_inout) {
#ifndef NDEBUG
    for (int32_t i=0; i<N_Values.N; i++) assert(fabs(real_inout[2*i])+fabs(imag_inout[2*i])<1e-20);
    for (int32_t i=0; i<N_Values.Ns2; i++) assert(fabs(imag_inout[2*i+1]+imag_inout[N_Values._2N-1-2*i])<1e-20);
#endif
}

EXPORT void IntPolynomial_ifft(LagrangeHalfCPolynomial* result, const IntPolynomial* p) {
    // fp1024_nayuki.execute_reverse_int(result->coefsC, p->coefs);
    // cplx* res          =>   result->coefsC
    // const int32_t* a   =>   p->coefs

    cplx* res = (cplx*)result->coefsC;
    const int32_t* a = p->coefs;

    double real_inout[N_Values._2N];
    double imag_inout[N_Values._2N];

    double* res_dbl=(double*) res;
    for (int32_t i=0; i<N_Values.N; i++) real_inout[i]=a[i]/2.;
    for (int32_t i=0; i<N_Values.N; i++) real_inout[N_Values.N+i]=-real_inout[i];
    for (int32_t i=0; i<N_Values._2N; i++) imag_inout[i]=0;
    check_alternate_real(real_inout, imag_inout);

    fft_transform_reverse(real_inout, imag_inout);

    for (int32_t i=0; i<N_Values.N; i+=2) {
        res_dbl[i]=real_inout[i+1];
        res_dbl[i+1]=imag_inout[i+1];
    }
    for (int32_t i=0; i<N_Values.Ns2; i++) {
	    assert(abs(cplx(real_inout[2*i+1],imag_inout[2*i+1])-res[i])<1e-20);
    }
    check_conjugate_cplx(real_inout, imag_inout);
}

EXPORT void TorusPolynomial_ifft(LagrangeHalfCPolynomial* result, const TorusPolynomial* p) {
    cplx* res = result->coefsC;
    Torus32* a = p->coefsT;

    double real_inout[N_Values._2N];
    double imag_inout[N_Values._2N];

    static const double _2pm33 = 1./double(INT64_C(1)<<33);
    int32_t* aa = (int32_t*) a;
    for (int32_t i=0; i<N_Values.N; i++) real_inout[i]=aa[i]*_2pm33;
    for (int32_t i=0; i<N_Values.N; i++) real_inout[N_Values.N+i]=-real_inout[i];
    for (int32_t i=0; i<N_Values._2N; i++) imag_inout[i]=0;
    check_alternate_real(real_inout, imag_inout);

    fft_transform_reverse(real_inout, imag_inout);

    for (int32_t i=0; i<N_Values.Ns2; i++) res[i]=cplx(real_inout[2*i+1],imag_inout[2*i+1]);
    check_conjugate_cplx(real_inout, imag_inout);
}

EXPORT void TorusPolynomial_fft(TorusPolynomial* result, const LagrangeHalfCPolynomial* p) {
    Torus32* res = result->coefsT;
    cplx* a = p->coefsC;

    double real_inout[N_Values._2N];
    double imag_inout[N_Values._2N];

    static const double _2p32 = double(INT64_C(1)<<32);
    static const double _1sN = double(1)/double(N_Values.N);
    //double* a_dbl=(double*) a;
    for (int32_t i=0; i<N_Values.N; i++) real_inout[2*i]=0;
    for (int32_t i=0; i<N_Values.N; i++) imag_inout[2*i]=0;
    for (int32_t i=0; i<N_Values.Ns2; i++) real_inout[2*i+1]=real(a[i]);
    for (int32_t i=0; i<N_Values.Ns2; i++) imag_inout[2*i+1]=imag(a[i]);
    for (int32_t i=0; i<N_Values.Ns2; i++) real_inout[N_Values._2N-1-2*i]=real(a[i]);
    for (int32_t i=0; i<N_Values.Ns2; i++) imag_inout[N_Values._2N-1-2*i]=-imag(a[i]);
#ifndef NDEBUG
    for (int32_t i=0; i<N_Values.N; i++) assert(real_inout[2*i]==0);
    for (int32_t i=0; i<N_Values.N; i++) assert(imag_inout[2*i]==0);
    for (int32_t i=0; i<N_Values.Ns2; i++) assert(real_inout[2*i+1]==real(a[i]));
    for (int32_t i=0; i<N_Values.Ns2; i++) assert(imag_inout[2*i+1]==imag(a[i]));
    for (int32_t i=0; i<N_Values.Ns2; i++) assert(real_inout[N_Values._2N-1-2*i]==real(a[i]));
    for (int32_t i=0; i<N_Values.Ns2; i++) assert(imag_inout[N_Values._2N-1-2*i]==-imag(a[i]));
    check_conjugate_cplx(real_inout, imag_inout);
#endif
    fft_transform(real_inout, imag_inout);
    for (int32_t i=0; i<N_Values.N; i++) res[i]=Torus32(int64_t(real_inout[i]*_1sN*_2p32));
    //pas besoin du fmod... Torus32(int64_t(fmod(rev_out[i]*_1sN,1.)*_2p32));
    check_alternate_real(real_inout, imag_inout);
}
