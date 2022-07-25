#include <cassert>
#include <cmath>
#include <cstdlib>
#include "tfhe_core.h"
#include "polynomials_arithmetic.h"
#include "lagrangehalfc_arithmetic.h"
#include "numeric_functions.h"

#include "fft.h"

using namespace std;

LagrangeHalfCPolynomial::LagrangeHalfCPolynomial(const int32_t N) {
    assert(N==1024);
    coefsC = new cplx[N/2];
    proc = &fp1024_nayuki;
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
    const int32_t Ns2 = reps->proc->Ns2;
    for (int32_t i=0; i<Ns2; i++)
	reps->coefsC[i] = 0;
}

EXPORT void LagrangeHalfCPolynomialSetTorusConstant(LagrangeHalfCPolynomial* result, const Torus32 mu) {
    const int32_t Ns2 = result->proc->Ns2;
    cplx* b = result->coefsC;
    const cplx muc = t32tod(mu);
    for (int32_t j=0; j<Ns2; j++)
    	b[j]=muc;
}

EXPORT void LagrangeHalfCPolynomialAddTorusConstant(LagrangeHalfCPolynomial* result, const Torus32 mu) {
    const int32_t Ns2 = result->proc->Ns2;
    cplx* b = result->coefsC;
    const cplx muc = t32tod(mu);
    for (int32_t j=0; j<Ns2; j++)
    	b[j]+=muc;
}

EXPORT void LagrangeHalfCPolynomialSetXaiMinusOne(LagrangeHalfCPolynomial* result, const int32_t ai) {
    const int32_t Ns2 = result->proc->Ns2;
    const int32_t _2N = result->proc->_2N;
    const cplx* omegaxminus1 = result->proc->omegaxminus1;
    for (int32_t i=0; i<Ns2; i++)
	result->coefsC[i]=omegaxminus1[((2*i+1)*ai)%_2N];
}

/** termwise multiplication in Lagrange space */
EXPORT void LagrangeHalfCPolynomialMul(
	LagrangeHalfCPolynomial* result,
	const LagrangeHalfCPolynomial* a,
	const LagrangeHalfCPolynomial* b) {
    const int32_t Ns2 = result->proc->Ns2;
    cplx* aa = a->coefsC;
    cplx* bb = b->coefsC;
    cplx* rr = result->coefsC;
    for (int32_t i=0; i<Ns2; i++)
	rr[i] = aa[i]*bb[i];
}

/** termwise multiplication and addTo in Lagrange space */
EXPORT void LagrangeHalfCPolynomialAddMul(
	LagrangeHalfCPolynomial* accum,
	const LagrangeHalfCPolynomial* a,
	const LagrangeHalfCPolynomial* b)
{
    const int32_t Ns2 = accum->proc->Ns2;
    cplx* aa = a->coefsC;
    cplx* bb = b->coefsC;
    cplx* rr = accum->coefsC;
    for (int32_t i=0; i<Ns2; i++)
	rr[i] += aa[i]*bb[i];
}

/** termwise multiplication and addTo in Lagrange space */
EXPORT void LagrangeHalfCPolynomialSubMul(
	LagrangeHalfCPolynomial* accum,
	const LagrangeHalfCPolynomial* a,
	const LagrangeHalfCPolynomial* b)
{
    const int32_t Ns2 = accum->proc->Ns2;
    cplx* aa = a->coefsC;
    cplx* bb = b->coefsC;
    cplx* rr = accum->coefsC;
    for (int32_t i=0; i<Ns2; i++)
	rr[i] -= aa[i]*bb[i];
}

EXPORT void LagrangeHalfCPolynomialAddTo(
	LagrangeHalfCPolynomial* accum,
	const LagrangeHalfCPolynomial* a) {
    const int32_t Ns2 = accum->proc->Ns2;
    cplx* aa = a->coefsC;
    cplx* rr = accum->coefsC;
    for (int32_t i=0; i<Ns2; i++)
	rr[i] += aa[i];
}

FFT_Processor_nayuki fp1024_nayuki(1024);

EXPORT void IntPolynomial_ifft(LagrangeHalfCPolynomial* result, const IntPolynomial* p) {
    fp1024_nayuki.execute_reverse_int(result->coefsC, p->coefs);
}
EXPORT void TorusPolynomial_ifft(LagrangeHalfCPolynomial* result, const TorusPolynomial* p) {
    fp1024_nayuki.execute_reverse_torus32(result->coefsC, p->coefsT);
}
EXPORT void TorusPolynomial_fft(TorusPolynomial* result, const LagrangeHalfCPolynomial* p) {
    fp1024_nayuki.execute_direct_torus32(result->coefsT, p->coefsC);
}

FFT_Processor_nayuki::FFT_Processor_nayuki(const int32_t N): _2N(2*N),N(N),Ns2(N/2) {
    real_inout = (double*) malloc(sizeof(double) * _2N);
    imag_inout = (double*) malloc(sizeof(double) * _2N);
    tables_direct = fft_init(_2N);
    tables_reverse = fft_init_reverse(_2N);
    omegaxminus1 = (cplx*) malloc(sizeof(cplx) * _2N);
    for (int32_t x=0; x<_2N; x++) {
        omegaxminus1[x]=cplx(cos(x*M_PI/N)-1., sin(x*M_PI/N)); // instead of cos(x*M_PI/N)-1. + sin(x*M_PI/N) * 1i
        //exp(i.x.pi/N)-1
    }


}

void FFT_Processor_nayuki::check_alternate_real() {
#ifndef NDEBUG
    for (int32_t i=0; i<_2N; i++) assert(fabs(imag_inout[i])<1e-8);
    for (int32_t i=0; i<N; i++) assert(fabs(real_inout[i]+real_inout[N+i])<1e-9);
#endif
}

void FFT_Processor_nayuki::check_conjugate_cplx() {
#ifndef NDEBUG
    for (int32_t i=0; i<N; i++) assert(fabs(real_inout[2*i])+fabs(imag_inout[2*i])<1e-20);
    for (int32_t i=0; i<Ns2; i++) assert(fabs(imag_inout[2*i+1]+imag_inout[_2N-1-2*i])<1e-20);
#endif
}

void FFT_Processor_nayuki::execute_reverse_int(cplx* res, const int32_t* a) {
    double* res_dbl=(double*) res;
    for (int32_t i=0; i<N; i++) real_inout[i]=a[i]/2.;
    for (int32_t i=0; i<N; i++) real_inout[N+i]=-real_inout[i];
    for (int32_t i=0; i<_2N; i++) imag_inout[i]=0;
    check_alternate_real();

    // fft_transform_reverse(tables_reverse,real_inout,imag_inout);
    fpga_fft_transform_reverse(tables_reverse, real_inout, imag_inout);

    for (int32_t i=0; i<N; i+=2) {
	res_dbl[i]=real_inout[i+1];
	res_dbl[i+1]=imag_inout[i+1];
    }
    for (int32_t i=0; i<Ns2; i++) {
	assert(abs(cplx(real_inout[2*i+1],imag_inout[2*i+1])-res[i])<1e-20);
    }
    check_conjugate_cplx();
}

void FFT_Processor_nayuki::execute_reverse_torus32(cplx* res, const Torus32* a) {
    static const double _2pm33 = 1./double(INT64_C(1)<<33);
    int32_t* aa = (int32_t*) a;
    for (int32_t i=0; i<N; i++) real_inout[i]=aa[i]*_2pm33;
    for (int32_t i=0; i<N; i++) real_inout[N+i]=-real_inout[i];
    for (int32_t i=0; i<_2N; i++) imag_inout[i]=0;
    check_alternate_real();

    // fft_transform_reverse(tables_reverse,real_inout,imag_inout);
    fpga_fft_transform_reverse(tables_reverse, real_inout, imag_inout);

    for (int32_t i=0; i<Ns2; i++) res[i]=cplx(real_inout[2*i+1],imag_inout[2*i+1]);
    check_conjugate_cplx();
}

void FFT_Processor_nayuki::execute_direct_torus32(Torus32* res, const cplx* a) {
    static const double _2p32 = double(INT64_C(1)<<32);
    static const double _1sN = double(1)/double(N);
    //double* a_dbl=(double*) a;
    for (int32_t i=0; i<N; i++) real_inout[2*i]=0;
    for (int32_t i=0; i<N; i++) imag_inout[2*i]=0;
    for (int32_t i=0; i<Ns2; i++) real_inout[2*i+1]=real(a[i]);
    for (int32_t i=0; i<Ns2; i++) imag_inout[2*i+1]=imag(a[i]);
    for (int32_t i=0; i<Ns2; i++) real_inout[_2N-1-2*i]=real(a[i]);
    for (int32_t i=0; i<Ns2; i++) imag_inout[_2N-1-2*i]=-imag(a[i]);
#ifndef NDEBUG
    for (int32_t i=0; i<N; i++) assert(real_inout[2*i]==0);
    for (int32_t i=0; i<N; i++) assert(imag_inout[2*i]==0);
    for (int32_t i=0; i<Ns2; i++) assert(real_inout[2*i+1]==real(a[i]));
    for (int32_t i=0; i<Ns2; i++) assert(imag_inout[2*i+1]==imag(a[i]));
    for (int32_t i=0; i<Ns2; i++) assert(real_inout[_2N-1-2*i]==real(a[i]));
    for (int32_t i=0; i<Ns2; i++) assert(imag_inout[_2N-1-2*i]==-imag(a[i]));
    check_conjugate_cplx();
#endif
    fft_transform(tables_direct,real_inout,imag_inout);
    for (int32_t i=0; i<N; i++) res[i]=Torus32(int64_t(real_inout[i]*_1sN*_2p32));
    //pas besoin du fmod... Torus32(int64_t(fmod(rev_out[i]*_1sN,1.)*_2p32));
    check_alternate_real();
}

FFT_Processor_nayuki::~FFT_Processor_nayuki() {
    fft_destroy(tables_direct);
    fft_destroy(tables_reverse);
    free(real_inout);
    free(imag_inout);
    free(omegaxminus1);
}

