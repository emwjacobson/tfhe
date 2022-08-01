#include <cassert>
#include <cmath>
#include <cstdlib>
#include "tfhe_core.h"
#include "polynomials_arithmetic.h"
#include "lagrangehalfc_arithmetic.h"
#include "numeric_functions.h"
#include "polynomials.h"
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
    LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3);
    IntPolynomial_ifft(tmp[0].coefsC, poly1->coefs);
    TorusPolynomial_ifft(tmp[1].coefsC, poly2->coefsT);
    LagrangeHalfCPolynomialMul(&tmp[2], &tmp[0], &tmp[1]);
    TorusPolynomial_fft(result->coefsT, tmp[2].coefsC);
    delete_LagrangeHalfCPolynomial_array(3,tmp);
}
EXPORT void torusPolynomialAddMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3);
    TorusPolynomial* tmpr = new_TorusPolynomial();
    IntPolynomial_ifft(tmp[0].coefsC,poly1->coefs);
    TorusPolynomial_ifft(tmp[1].coefsC, poly2->coefsT);
    LagrangeHalfCPolynomialMul(&tmp[2], &tmp[0], &tmp[1]);
    TorusPolynomial_fft(tmpr->coefsT, tmp[2].coefsC);
    torusPolynomialAddTo(result, tmpr);
    delete_TorusPolynomial(tmpr);
    delete_LagrangeHalfCPolynomial_array(3,tmp);
}
EXPORT void torusPolynomialSubMulRFFT(TorusPolynomial* result, const IntPolynomial* poly1, const TorusPolynomial* poly2) {
    LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3);
    TorusPolynomial* tmpr = new_TorusPolynomial();
    IntPolynomial_ifft(tmp[0].coefsC, poly1->coefs);
    TorusPolynomial_ifft(tmp[1].coefsC, poly2->coefsT);
    LagrangeHalfCPolynomialMul(&tmp[2], &tmp[0], &tmp[1]);
    TorusPolynomial_fft(tmpr->coefsT, tmp[2].coefsC);
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

EXPORT void IntPolynomial_ifft(LagrangeHalfCPolynomial_Collapsed result, const IntPolynomial_Collapsed p) {
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

	memcpy(result, result_map, sizeof(cplx) * Value_Ns2);
}

EXPORT void TorusPolynomial_ifft(LagrangeHalfCPolynomial_Collapsed result, const TorusPolynomial_Collapsed p) {
    cl::Buffer result_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(LagrangeHalfCPolynomial_Collapsed));
    cl::Buffer p_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TorusPolynomial_Collapsed));

	cplx *result_map = (cplx *)fpga.q.enqueueMapBuffer(result_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(LagrangeHalfCPolynomial_Collapsed));
	Torus32 *p_map = (Torus32 *)fpga.q.enqueueMapBuffer(p_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TorusPolynomial_Collapsed));

	memcpy(p_map, p, sizeof(int32_t) * Value_N);

	fpga.k_TorusPolynomial_ifft.setArg(0, result_buf);
	fpga.k_TorusPolynomial_ifft.setArg(1, p_buf);

	fpga.q.enqueueMigrateMemObjects({ result_buf, p_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_TorusPolynomial_ifft);
	fpga.q.enqueueMigrateMemObjects({ result_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

	memcpy(result, result_map, sizeof(cplx) * Value_Ns2);
}

EXPORT void TorusPolynomial_fft(TorusPolynomial_Collapsed result, const LagrangeHalfCPolynomial_Collapsed p) {
    cl::Buffer result_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(TorusPolynomial_Collapsed));
    cl::Buffer p_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(LagrangeHalfCPolynomial_Collapsed));

	Torus32 *result_map = (Torus32 *)fpga.q.enqueueMapBuffer(result_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(TorusPolynomial_Collapsed));
	cplx *p_map = (cplx *)fpga.q.enqueueMapBuffer(p_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(LagrangeHalfCPolynomial_Collapsed));

	memcpy(p_map, p, sizeof(cplx) * Value_Ns2);

	fpga.k_TorusPolynomial_fft.setArg(0, result_buf);
	fpga.k_TorusPolynomial_fft.setArg(1, p_buf);

	fpga.q.enqueueMigrateMemObjects({ result_buf, p_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_TorusPolynomial_fft);
	fpga.q.enqueueMigrateMemObjects({ result_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

	memcpy(result, result_map, sizeof(Torus32) * Value_N);
}
