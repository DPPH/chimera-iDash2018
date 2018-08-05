#ifndef FHE_BIGCOMPLEX_H
#define FHE_BIGCOMPLEX_H

#include <mpfr.h>
#include "BigReal.h"

/**
 * BigComplex represent fixed precision complex numbers of module <=1 having a no-overflow guarantee.
 * (for all arithmetic operations they are involved in, the module stays <=1)
 *
 * A BigComplex x of precision 2^-64.n is represented by two BigReal of same precision
 */

class BigComplex {
public:
    BigReal real;
    BigReal imag;

    BigComplex(uint64_t nblimbs);

    ~BigComplex();
};

/**
 * reference to a BigComplex
 */
class BigComplexRef {
public:
    BigReal *const real;
    BigReal *const imag;

    BigComplexRef(BigReal *real, BigReal *imag);

    BigComplexRef(const BigReal *real, const BigReal *imag);

    BigComplexRef(BigComplex &z);

    BigComplexRef(const BigComplex &z);
};

void accurate_power_unity(BigComplexRef dest, int i, int n);

void mul(BigComplexRef dest, BigComplexRef a, BigComplexRef b);

void mulTo(BigComplexRef dest, BigComplexRef a);

void add(BigComplexRef dest, BigComplexRef a, BigComplexRef b);

void sub(BigComplexRef dest, BigComplexRef a, BigComplexRef b);

void copy(BigComplexRef dest, BigComplexRef a);

BigComplex *new_BigComplex_array(uint64_t n, uint64_t nblimbs);

void delete_BigComplex_array(uint64_t n, BigComplex *array);


#endif //FHE_BIGCOMPLEX_H
