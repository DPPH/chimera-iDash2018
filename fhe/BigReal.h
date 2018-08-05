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

/** addition (warning, you must ensure that the result does not overflows) */
void add(BigReal &dest, const BigReal &a, const BigReal &b);

/** subtraction (warning, you must ensure that the result does not overflows) */
void sub(BigReal &dest, const BigReal &a, const BigReal &b);

/** negation */
void neg(BigReal &dest, const BigReal &a);

/** multiplication with an integer (warning, you must ensure that the result does not overflows) */
void extmul(BigReal &dest, int a, const BigReal &b);

/** multiplication */
void mul(BigReal &dest, const BigReal &a, const BigReal &b);

/** division by power of 2 */
void div2ui(BigReal &dest, const BigReal &a, uint64_t b);

/** copy */
void copy(BigReal &dest, const BigReal &a);

/** conversion torus to bigreal */
void to_BigReal(BigReal &dest, const BigTorusRef &v);

/** conversion bigreal to RR */
NTL::RR to_RR(const BigReal &v);

/** conversion RR to bigreal */
void to_BigReal(BigReal &dest, const NTL::RR &v);

/** create an array of BigReal */
BigReal *new_BigReal_array(uint64_t n, uint64_t nblimbs);

/** delete an array of BigReal */
void delete_BigReal_array(uint64_t n, BigReal *array);

#endif //FHE_BIG_REAL_H
