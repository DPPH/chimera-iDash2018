#ifndef FHE_BIG_REAL_H
#define FHE_BIG_REAL_H


#include <cstdint>
#include <gmp.h>
#include "commons.h"
#include "BigTorus.h"

/**
 * BigReal represent fixed precision real numbers between [-1,1] having a no-overflow guarantee.
 * (for all arithmetic operations they are involved in, their value stay between [-1,1])
 *
 * A BigReal x of precision 2^-64.n is represented by a the GMP integer 2^64.x
 */
class BigReal {
public:
    mpz_t value;
    uint64_t nblimbs;

    BigReal(uint64_t nblimbs);

    ~BigReal();

    BigReal(const BigReal &) = delete;

    void operator=(const BigReal &)= delete;
};

void add(BigReal &dest, const BigReal &a, const BigReal &b);

void sub(BigReal &dest, const BigReal &a, const BigReal &b);

void neg(BigReal &dest, const BigReal &a);

void extmul(BigReal &dest, int a, const BigReal &b);

void mul(BigReal &dest, const BigReal &a, const BigReal &b);

void div2ui(BigReal &dest, const BigReal &a, uint64_t b);

void copy(BigReal &dest, const BigReal &a);

void to_BigReal(BigReal &dest, const BigTorusRef &v);

//Real96 dtor96(double d) {
//    Real96 reps;
//    mpfr_set_d(reps.v, d, MPFR_RNDN);
//    return reps;
//}


#endif //FHE_BIG_REAL_H
