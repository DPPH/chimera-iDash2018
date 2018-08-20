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
    BigReal real; ///< real part
    BigReal imag; ///< imaginary part

    /* constructs a BigComplex of precision nblimbs x 64 */
    BigComplex(UINT64 nblimbs);

    /* destroys a BigComplex */
    ~BigComplex();

    /* prevents auto-copying */
    BigComplex(const BigComplex &) = delete;

    /* prevents auto-copying */
    void operator=(const BigComplex &)= delete;
};

/**
 * reference to a BigComplex (this only holds two pointers to existing real and imag)
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


/** serialize:
 *  magic number:   int64 on 8 bytes
 *  value:   2*BigReal
 *
 */
void serializeBigcomplexContent(std::ostream &out, const BigComplexRef &value);


/** serialize:
 *  magic number:   int64 on 8 bytes
 *  value:   2*BigReal
 *
 */
void deserializeBigComplexContent(std::istream &in, BigComplexRef &value);


/** serialize:
 * value:   2*BigReal
 */
void serializeBigComplex(std::ostream &out, const BigComplex &value);


/** serialize:
 * value:   2*BigReal
 *
 */
std::shared_ptr<BigComplex> deserializeBigComplex(std::istream &in);


/** computes exp(2i.pi. k/n) */
void accurate_power_unity(BigComplexRef dest, int k, int n);

/** multiplication dest = a * b */
void mul(BigComplexRef dest, BigComplexRef a, BigComplexRef b);

/** multiplication dest *= a */
void mulTo(BigComplexRef dest, BigComplexRef a);

/** addition (warning, you must ensure that this does not cause an overflow) */
void add(BigComplexRef dest, BigComplexRef a, BigComplexRef b);

/** subtraction (warning, you must ensure that this does not cause an overflow) */
void sub(BigComplexRef dest, BigComplexRef a, BigComplexRef b);

/** copy dest = a */
void copy(BigComplexRef dest, BigComplexRef a);

/** dest = 0 */
void zero(BigComplexRef dest);

/** dest = dest + a * b */
void addMulTo(BigComplexRef dest, const BigComplexRef &a, const BigComplexRef &b);

/** create an array */
BigComplex *new_BigComplex_array(UINT64 n, UINT64 nblimbs);

/** delete an array */
void delete_BigComplex_array(UINT64 n, BigComplex *array);


#endif //FHE_BIGCOMPLEX_H
