#include <gmp.h>
#include "BigTorus.h"
#include "commons.h"
#include "arithmetic.h"

BigTorusParams::BigTorusParams(uint64_t torus_limbs) : torus_limbs(torus_limbs) {}

BigTorusRef::BigTorusRef(uint64_t *limbs_raw, const BigTorusParams *params) : limbs_raw(limbs_raw), params(params) {}

BigTorusRef::BigTorusRef(const BigTorus &torus) : limbs_raw(torus.limbs_raw), params(torus.params) {}

BigTorusRef::BigTorusRef(BigTorus &torus) : limbs_raw(torus.limbs_raw), params(torus.params) {}

void bigTorusRawScale(uint64_t *limbs, int64_t coef, uint64_t nblimbs) {
    assert_dramatically(coef > 0, "negative scale not supported");
    mpn_mul_1(limbs, limbs, nblimbs, coef);
}

void bigTorusScale(const BigTorusRef &x, int64_t coef) {
    bigTorusRawScale(x.limbs_raw, coef, x.params->torus_limbs);
}

void copy(BigTorusRef dest, const BigTorusRef &x, uint64_t limb_precision) {
    uint64_t xsize = x.params->torus_limbs;
    uint64_t dsize = dest.params->torus_limbs;
    if (limb_precision == NA) {
        limb_precision = std::min(xsize, dsize);
    } else {
        assert_dramatically(limb_precision <= dsize, "destination is not precise enough");
        assert_dramatically(limb_precision <= xsize, "source is not precise enough");
    }
    uint64_t xoffset = xsize - limb_precision;
    uint64_t doffset = dsize - limb_precision;
    //copy the significant bits
    mpn_copyi(dest.limbs_raw + doffset, x.limbs_raw + xoffset, limb_precision);
}

void random(BigTorusRef dest, uint64_t limb_precision) {
    uint64_t dsize = dest.params->torus_limbs;
    if (limb_precision == NA) {
        limb_precision = dsize;
    } else {
        assert_dramatically(limb_precision <= dsize, "destination is not precise enough");
    }
    uint64_t doffset = dsize - limb_precision;
    mpn_random(dest.limbs_raw + doffset, limb_precision);
}

void add(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, uint64_t limb_precision) {
    uint64_t asize = a.params->torus_limbs;
    uint64_t bsize = b.params->torus_limbs;
    uint64_t dsize = dest.params->torus_limbs;
    if (limb_precision == NA) {
        limb_precision = std::min(std::min(asize, bsize), dsize);
    } else {
        assert_dramatically(limb_precision <= dsize, "destination is not precise enough");
        assert_dramatically(limb_precision <= asize, "source is not precise enough");
        assert_dramatically(limb_precision <= bsize, "source is not precise enough");
    }
    uint64_t aoffset = asize - limb_precision;
    uint64_t boffset = asize - limb_precision;
    uint64_t doffset = dsize - limb_precision;
    mpn_add_n(dest.limbs_raw + doffset, a.limbs_raw + aoffset, b.limbs_raw + boffset, limb_precision);
}

void sub(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, uint64_t limb_precision) {
    uint64_t asize = a.params->torus_limbs;
    uint64_t bsize = b.params->torus_limbs;
    uint64_t dsize = dest.params->torus_limbs;
    if (limb_precision == NA) {
        limb_precision = std::min(std::min(asize, bsize), dsize);
    } else {
        assert_dramatically(limb_precision <= dsize, "destination is not precise enough");
        assert_dramatically(limb_precision <= asize, "source is not precise enough");
        assert_dramatically(limb_precision <= bsize, "source is not precise enough");
    }
    uint64_t aoffset = asize - limb_precision;
    uint64_t boffset = asize - limb_precision;
    uint64_t doffset = dsize - limb_precision;
    mpn_sub_n(dest.limbs_raw + doffset, a.limbs_raw + aoffset, b.limbs_raw + boffset, limb_precision);
}

void add_noise(BigTorusRef dest, uint64_t alpha_bits, uint64_t limb_precision) {
    uint64_t dsize = dest.params->torus_limbs;
    if (limb_precision == NA) {
        limb_precision = dsize;
    } else {
        assert_dramatically(limb_precision <= dsize, "destination is not precise enough");
    }
    //uint64_t doffset = dsize - limb_precision;
    uint64_t i;
    uint64_t ab;
    for (i = 0, ab = alpha_bits; i < limb_precision; i++, ab -= BITS_PER_LIMBS) {
        if (ab >= BITS_PER_LIMBS) continue; //nothing
        uint64_t r = random_uint64_t();
        if (ab > 0) { //between 1 and 63
            uint64_t mask = 1ul << (64ul - ab);
            mask |= (mask - 1);
            r &= mask;
        }
        dest.limbs_raw[dsize - 1 - i] ^= r;
    }
}

