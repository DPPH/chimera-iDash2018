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
    UINT64 nblimbs;

    BigReal(UINT64 nblimbs);

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
void div2ui(BigReal &dest, const BigReal &a, UINT64 b);

/** copy */
void copy(BigReal &dest, const BigReal &a);

/** zero */
void zero(BigReal &dest);

/** conversion torus to bigreal (outputs v between -1/2 and 1/2) */
void to_BigReal(BigReal &dest, const BigTorusRef &v);

/** conversion torus to bigreal (outputs a/2^nbits) */
void to_BigReal(BigReal &dest, const int64_t a, UINT64 a_nbits);

/** conversion bigreal to torus (outputs a*2^lshiftbits mod 1 with given limb precision) */
void to_BigTorus(BigTorusRef dest, const BigReal &a, UINT64 lshift_bits, UINT64 out_precision_limbs);

/** conversion bigreal to RR */
NTL::RR to_RR(const BigReal &v);

/** conversion RR to bigreal */
void to_BigReal(BigReal &dest, const NTL::RR &v);

/** create an array of BigReal */
BigReal *new_BigReal_array(UINT64 n, UINT64 nblimbs);

/** delete an array of BigReal */
void delete_BigReal_array(UINT64 n, BigReal *array);

#endif //FHE_BIG_REAL_H
