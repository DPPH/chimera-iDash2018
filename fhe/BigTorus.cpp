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

bigtorus_addsub_params
prepare_addsub(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, uint64_t limb_precision) {
    bigtorus_addsub_params r;
    r.asize = a.params->torus_limbs;
    r.bsize = b.params->torus_limbs;
    r.dsize = dest.params->torus_limbs;
    if (limb_precision == NA) {
        r.limb_precision = std::min(std::min(r.asize, r.bsize), r.dsize);
    } else {
        assert_dramatically(limb_precision <= r.dsize, "destination is not precise enough");
        assert_dramatically(limb_precision <= r.asize, "source is not precise enough");
        assert_dramatically(limb_precision <= r.bsize, "source is not precise enough");
        r.limb_precision = limb_precision;
    }
    r.aoffset = r.asize - r.limb_precision;
    r.boffset = r.bsize - r.limb_precision;
    r.doffset = r.dsize - r.limb_precision;
    return r;
}

void add_prep(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, const bigtorus_addsub_params &prep) {
    mpn_add_n(dest.limbs_raw + prep.doffset,
              a.limbs_raw + prep.aoffset,
              b.limbs_raw + prep.boffset,
              prep.limb_precision);
}

void sub_prep(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, const bigtorus_addsub_params &prep) {
    mpn_sub_n(dest.limbs_raw + prep.doffset,
              a.limbs_raw + prep.aoffset,
              b.limbs_raw + prep.boffset,
              prep.limb_precision);
}

void add(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, uint64_t limb_precision) {
    add_prep(dest, a, b, prepare_addsub(dest, a, b, limb_precision));
}

void sub(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, uint64_t limb_precision) {
    sub_prep(dest, a, b, prepare_addsub(dest, a, b, limb_precision));
}

void zero(BigTorusRef dest) {
    mpn_zero(dest.limbs_raw, dest.params->torus_limbs);
}
