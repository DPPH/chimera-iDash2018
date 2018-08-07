#include <cassert>
#include "BigReal.h"

BigReal::BigReal(uint64_t nblimbs) : nblimbs(nblimbs) {
    mpz_init(value);
}

BigReal::~BigReal() {
    mpz_clear(value);
}

void add(BigReal &dest, const BigReal &a, const BigReal &b) {
    mpz_add(dest.value, a.value, b.value);
}

void sub(BigReal &dest, const BigReal &a, const BigReal &b) {
    mpz_sub(dest.value, a.value, b.value);
}

void neg(BigReal &dest, const BigReal &a) {
    mpz_neg(dest.value, a.value);
}

void extmul(BigReal &dest, int a, const BigReal &b) {
    mpz_mul_si(dest.value, b.value, a);
}

void mul(BigReal &dest, const BigReal &a, const BigReal &b) {
    mpz_mul(dest.value, a.value, b.value);
    mpz_cdiv_q_2exp(dest.value, dest.value, dest.nblimbs * BITS_PER_LIMBS);
}

void div2ui(BigReal &dest, const BigReal &a, uint64_t b) {
    mpz_cdiv_q_2exp(dest.value, a.value, b);
}

void to_BigReal(BigReal &dest, const BigTorusRef &v) {
    const int dsize = dest.nblimbs;
    const int vsize = v.params->torus_limbs;
    const int vcopysize = std::min(dsize, vsize);
    const int zerofill = dsize - vcopysize;
    const int vcopyoffset = vsize - vcopysize;
    //copy the nblimbs most significant limbs
    mpz_realloc2(dest.value, dsize * BITS_PER_LIMBS);
    assert(dest.value->_mp_alloc >= dsize);
    //dest.value->_mp_d = (mp_limb_t *) malloc(dsize * BYTES_PER_LIMBS);
    //fill the zeros
    if (zerofill > 0) {
        mpn_zero(dest.value->_mp_d, zerofill);
    }
    //negate or copy limbs from v
    bool vnegative = v.limbs[vsize - 1] >> 63u;
    if (vnegative) {
        //v is negative
        mpn_neg(dest.value->_mp_d + zerofill, v.limbs + vcopyoffset, vcopysize);
    } else {
        //v is positive
        mpn_copyi(dest.value->_mp_d + zerofill, v.limbs + vcopyoffset, vcopysize);
    }
    //find the right size
    int dl;
    for (dl = dsize; dl >= 1; dl--) {
        if (dest.value->_mp_d[dl - 1] != 0) break;
    }
    dest.value->_mp_size = vnegative ? -dl : dl;
}

void copy(BigReal &dest, const BigReal &a) {
    mpz_set(dest.value, a.value);
}

BigReal *new_BigReal_array(uint64_t n, uint64_t nblimbs) {
    BigReal *reps = (BigReal *) malloc(n * sizeof(BigReal));
    for (uint64_t i = 0; i < n; i++) {
        new(reps + i) BigReal(nblimbs);
    }
    return reps;
}

void delete_BigReal_array(uint64_t n, BigReal *array) {
    for (uint64_t i = 0; i < n; i++) {
        (array + i)->~BigReal();
    }
    free(array);
}

#include <NTL/ZZ_limbs.h>
#include <gmp.h>

using namespace NTL;

NTL::RR to_RR(const BigReal &v) {
    ZZ vv;
    if (v.value->_mp_size == 0) {
        clear(vv);
    } else if (v.value->_mp_size > 0) {
        ZZ_limbs_set(vv, v.value->_mp_d, v.value->_mp_size);
    } else {
        ZZ_limbs_set(vv, v.value->_mp_d, -v.value->_mp_size);
        negate(vv, vv);
    }
    return to_RR(vv) / pow(to_RR(2), to_RR(v.nblimbs * BITS_PER_LIMBS));
}

void to_BigReal(BigReal &dest, const NTL::RR &v) {
    assert_dramatically(abs(v) <= 1, "Bad range for BigReal");
    RR::SetPrecision(dest.nblimbs * BITS_PER_LIMBS);
    ZZ vv = RoundToZZ(v * pow(to_RR(2), to_RR(dest.nblimbs * BITS_PER_LIMBS)));
    bool vneg = vv < 0;
    int64_t vsize = vv.size();
    if (vneg) negate(vv, vv);
    if (vv == 0) {
        mpz_set_ui(dest.value, 0);
    } else {
        mpz_realloc2(dest.value, vsize * BITS_PER_LIMBS);
        assert(dest.value->_mp_alloc >= vsize);
        dest.value->_mp_size = vneg ? -vsize : vsize;
        //dest.value->_mp_d = (mp_limb_t *) malloc(vsize * BYTES_PER_LIMBS);
        mpn_copyi(dest.value->_mp_d, ZZ_limbs_get(vv), vsize);
    }
}
