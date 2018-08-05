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
    mpz_clear(dest.value);
    dest.value->_mp_alloc = dsize;
    dest.value->_mp_d = (mp_limb_t *) malloc(dsize * BYTES_PER_LIMBS);
    //fill the zeros
    if (zerofill > 0) {
        mpn_zero(dest.value->_mp_d, zerofill);
    }
    //negate or copy limbs from v
    bool vnegative = v.limbs_raw[vsize - 1] >> 63u;
    if (vnegative) {
        //v is negative
        mpn_neg(dest.value->_mp_d + zerofill, v.limbs_raw + vcopyoffset, vcopysize);
    } else {
        //v is positive
        mpn_copyi(dest.value->_mp_d + zerofill, v.limbs_raw + vcopyoffset, vcopysize);
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
