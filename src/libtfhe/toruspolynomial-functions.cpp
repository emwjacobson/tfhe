#include <cassert>
#include <cmath>
#include "tfhe_core.h"
#include "numeric_functions.h"
#include "polynomials.h"

using namespace std;

// TorusPolynomial = 0
EXPORT void torusPolynomialClear(TorusPolynomial *result) {
    for (int32_t i = 0; i < Value_N; ++i) result->coefsT[i] = 0;
}

// TorusPolynomial = random
EXPORT void torusPolynomialUniform(TorusPolynomial *result) {
    Torus32 *x = result->coefsT;

    for (int32_t i = 0; i < Value_N; ++i)
        x[i] = uniformTorus32_distrib(generator);
}

// TorusPolynomial = TorusPolynomial
EXPORT void torusPolynomialCopy(
        TorusPolynomial *result,
        const TorusPolynomial *sample) {
    assert(result != sample);
    const Torus32 *__restrict s = sample->coefsT;
    Torus32 *__restrict r = result->coefsT;

    for (int32_t i = 0; i < Value_N; ++i) r[i] = s[i];
}

// TorusPolynomial + TorusPolynomial
EXPORT void torusPolynomialAdd(TorusPolynomial *result, const TorusPolynomial *poly1, const TorusPolynomial *poly2) {
    assert(result != poly1); //if it fails here, please use addTo
    assert(result != poly2); //if it fails here, please use addTo
    Torus32 *__restrict r = result->coefsT;
    const Torus32 *__restrict a = poly1->coefsT;
    const Torus32 *__restrict b = poly2->coefsT;

    for (int32_t i = 0; i < Value_N; ++i)
        r[i] = a[i] + b[i];
}

// TorusPolynomial += TorusPolynomial
EXPORT void torusPolynomialAddTo(TorusPolynomial *result, const TorusPolynomial *poly2) {
    Torus32 *r = result->coefsT;
    const Torus32 *b = poly2->coefsT;

    for (int32_t i = 0; i < Value_N; ++i)
        r[i] += b[i];
}


// TorusPolynomial - TorusPolynomial
EXPORT void torusPolynomialSub(TorusPolynomial *result, const TorusPolynomial *poly1, const TorusPolynomial *poly2) {
    assert(result != poly1); //if it fails here, please use subTo
    assert(result != poly2); //if it fails here, please use subTo
    Torus32 *__restrict r = result->coefsT;
    const Torus32 *a = poly1->coefsT;
    const Torus32 *b = poly2->coefsT;

    for (int32_t i = 0; i < Value_N; ++i)
        r[i] = a[i] - b[i];
}

// TorusPolynomial -= TorusPolynomial
EXPORT void torusPolynomialSubTo(TorusPolynomial *result, const TorusPolynomial *poly2) {
    Torus32 *r = result->coefsT;
    const Torus32 *b = poly2->coefsT;

    for (int32_t i = 0; i < Value_N; ++i)
        r[i] -= b[i];
}

// TorusPolynomial + p*TorusPolynomial
EXPORT void
torusPolynomialAddMulZ(TorusPolynomial *result, const TorusPolynomial *poly1, int32_t p, const TorusPolynomial *poly2) {
    Torus32 *r = result->coefsT;
    const Torus32 *a = poly1->coefsT;
    const Torus32 *b = poly2->coefsT;

    for (int32_t i = 0; i < Value_N; ++i)
        r[i] = a[i] + p * b[i];
}

// TorusPolynomial += p*TorusPolynomial
EXPORT void torusPolynomialAddMulZTo(TorusPolynomial *result, const int32_t p, const TorusPolynomial *poly2) {
    Torus32 *r = result->coefsT;
    const Torus32 *b = poly2->coefsT;

    for (int32_t i = 0; i < Value_N; ++i) r[i] += p * b[i];
}

// TorusPolynomial - p*TorusPolynomial
EXPORT void torusPolynomialSubMulZ(TorusPolynomial *result, const TorusPolynomial *poly1, const int32_t p,
                                   const TorusPolynomial *poly2) {
    Torus32 *r = result->coefsT;
    const Torus32 *a = poly1->coefsT;
    const Torus32 *b = poly2->coefsT;

    for (int32_t i = 0; i < Value_N; ++i) r[i] = a[i] - p * b[i];
}

//result= (X^{a}-1)*source
EXPORT void torusPolynomialMulByXaiMinusOne(TorusPolynomial *result, int32_t a, const TorusPolynomial *source) {
    Torus32 *out = result->coefsT;
    const Torus32 *in = source->coefsT;

    assert(a >= 0 && a < 2 * Value_N);

    if (a < Value_N) {
        for (int32_t i = 0; i < a; i++)//sur que i-a<0
            out[i] = -in[i - a + Value_N] - in[i];
        for (int32_t i = a; i < Value_N; i++)//sur que N>i-a>=0
            out[i] = in[i - a] - in[i];
    } else {
        const int32_t aa = a - Value_N;
        for (int32_t i = 0; i < aa; i++)//sur que i-a<0
            out[i] = in[i - aa + Value_N] - in[i];
        for (int32_t i = aa; i < Value_N; i++)//sur que N>i-a>=0
            out[i] = -in[i - aa] - in[i];
    }
}


//result= X^{a}*source
EXPORT void torusPolynomialMulByXai(TorusPolynomial *result, int32_t a, const TorusPolynomial *source) {
    Torus32 *out = result->coefsT;
    const Torus32 *in = source->coefsT;

    assert(a >= 0 && a < 2 * Value_N);
    assert(result != source);

    if (a < Value_N) {
        for (int32_t i = 0; i < a; i++)//sur que i-a<0
            out[i] = -in[i - a + Value_N];
        for (int32_t i = a; i < Value_N; i++)//sur que N>i-a>=0
            out[i] = in[i - a];
    } else {
        const int32_t aa = a - Value_N;
        for (int32_t i = 0; i < aa; i++)//sur que i-a<0
            out[i] = in[i - aa + Value_N];
        for (int32_t i = aa; i < Value_N; i++)//sur que N>i-a>=0
            out[i] = -in[i - aa];
    }
}


// TorusPolynomial -= p*TorusPolynomial
EXPORT void torusPolynomialSubMulZTo(TorusPolynomial *result, int32_t p, const TorusPolynomial *poly2) {
    Torus32 *r = result->coefsT;
    const Torus32 *b = poly2->coefsT;

    for (int32_t i = 0; i < Value_N; ++i) r[i] -= p * b[i];
}


// Norme Euclidienne d'un IntPolynomial
EXPORT double intPolynomialNormSq2(const IntPolynomial *poly) {
    int32_t temp1 = 0;

    for (int32_t i = 0; i < Value_N; ++i) {
        int32_t temp0 = poly->coefs[i] * poly->coefs[i];
        temp1 += temp0;
    }
    return temp1;
}

// Sets to zero
EXPORT void intPolynomialClear(IntPolynomial *poly) {
    for (int32_t i = 0; i < Value_N; ++i)
        poly->coefs[i] = 0;
}

// Sets to zero
EXPORT void intPolynomialCopy(IntPolynomial *result, const IntPolynomial *source) {
    for (int32_t i = 0; i < Value_N; ++i)
        result->coefs[i] = source->coefs[i];
}

/** accum += source */
EXPORT void intPolynomialAddTo(IntPolynomial *accum, const IntPolynomial *source) {
    for (int32_t i = 0; i < Value_N; ++i)
        accum->coefs[i] += source->coefs[i];
}

/**  result = (X^ai-1) * source */
EXPORT void intPolynomialMulByXaiMinusOne(IntPolynomial *result, int32_t ai, const IntPolynomial *source) {
    int32_t *out = result->coefs;
    const int32_t *in = source->coefs;

    assert(ai >= 0 && ai < 2 * Value_N);

    if (ai < Value_N) {
        for (int32_t i = 0; i < ai; i++)//sur que i-a<0
            out[i] = -in[i - ai + Value_N] - in[i];
        for (int32_t i = ai; i < Value_N; i++)//sur que N>i-a>=0
            out[i] = in[i - ai] - in[i];
    } else {
        const int32_t aa = ai - Value_N;
        for (int32_t i = 0; i < aa; i++)//sur que i-a<0
            out[i] = in[i - aa + Value_N] - in[i];
        for (int32_t i = aa; i < Value_N; i++)//sur que N>i-a>=0
            out[i] = -in[i - aa] - in[i];
    }
}



// Norme infini de la distance entre deux TorusPolynomial
EXPORT double torusPolynomialNormInftyDist(const TorusPolynomial *poly1, const TorusPolynomial *poly2) {
    double norm = 0;

    // Max between the coefficients of abs(poly1-poly2)
    for (int32_t i = 0; i < Value_N; ++i) {
        double r = abs(t32tod(poly1->coefsT[i] - poly2->coefsT[i]));
        if (r > norm) { norm = r; }
    }
    return norm;
}






// Norme 2 d'un IntPolynomial
EXPORT double intPolynomialNorm2sq(const IntPolynomial *poly) {
    double norm = 0;

    for (int32_t i = 0; i < Value_N; ++i) {
        double r = poly->coefs[i];
        norm += r * r;
    }
    return norm;
}

// Norme infini de la distance entre deux IntPolynomial
EXPORT double intPolynomialNormInftyDist(const IntPolynomial *poly1, const IntPolynomial *poly2) {
    double norm = 0;


    // Max between the coefficients of abs(poly1-poly2)
    for (int32_t i = 0; i < Value_N; ++i) {
        double r = abs(poly1->coefs[i] - poly2->coefs[i]);
        if (r > norm) { norm = r; }
    }
    return norm;
}


