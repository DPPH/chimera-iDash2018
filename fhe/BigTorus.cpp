#include <gmp.h>
#include "BigTorus.h"
#include "commons.h"
#include "arithmetic.h"

BigTorusParams::BigTorusParams(UINT64 torus_limbs, int64_t plaintext_expo, int64_t level_expo) :
        torus_limbs(torus_limbs),
        plaintext_expo(plaintext_expo),
        level_expo(level_expo) {
}

BigTorusRef::BigTorusRef(UINT64 *limbs, const BigTorusParams *params) : limbs(limbs), params(params) {}

BigTorusRef::BigTorusRef(const BigTorus &torus) : limbs(torus.limbs), params(torus.params) {}

BigTorusRef::BigTorusRef(BigTorus &torus) : limbs(torus.limbs), params(torus.params) {}

void bigTorusRawScale(UINT64 *limbs, int64_t coef, UINT64 nblimbs) {
    assert_dramatically(coef > 0, "negative scale not supported");
    mpn_mul_1(limbs, limbs, nblimbs, coef);
}

void bigTorusScale(const BigTorusRef &x, int64_t coef) {
    bigTorusRawScale(x.limbs, coef, x.params->torus_limbs);
}

void copy(BigTorusRef dest, const BigTorusRef &x, UINT64 limb_precision) {
    UINT64 xsize = x.params->torus_limbs;
    UINT64 dsize = dest.params->torus_limbs;
    if (limb_precision == NA) {
        limb_precision = std::min(xsize, dsize);
    } else {
        assert_dramatically(limb_precision <= dsize, "destination is not precise enough");
        assert_dramatically(limb_precision <= xsize, "source is not precise enough");
    }
    UINT64 xoffset = xsize - limb_precision;
    UINT64 doffset = dsize - limb_precision;
    //copy the significant bits
    mpn_copyi(dest.limbs + doffset, x.limbs + xoffset, limb_precision);
}

void random(BigTorusRef dest, UINT64 limb_precision) {
    UINT64 dsize = dest.params->torus_limbs;
    if (limb_precision == NA) {
        limb_precision = dsize;
    } else {
        assert_dramatically(limb_precision <= dsize, "destination is not precise enough");
    }
    UINT64 doffset = dsize - limb_precision;
    mpn_random(dest.limbs + doffset, limb_precision);
}

void add_noise(BigTorusRef dest, UINT64 alpha_bits, UINT64 limb_precision) {
    UINT64 dsize = dest.params->torus_limbs;
    if (limb_precision == NA) {
        limb_precision = dsize;
    } else {
        assert_dramatically(limb_precision <= dsize, "destination is not precise enough");
    }
    //UINT64 doffset = dsize - limb_precision;
    UINT64 i;
    UINT64 ab;
    for (i = 0, ab = alpha_bits; i < limb_precision; i++, ab -= BITS_PER_LIMBS) {
        if (ab >= BITS_PER_LIMBS) continue; //nothing
        UINT64 r = random_uint64_t();
        if (ab > 0) { //between 1 and 63
            UINT64 mask = 1ul << (64ul - ab);
            mask |= (mask - 1);
            r &= mask;
        }
        dest.limbs[dsize - 1 - i] ^= r;
    }
}

bigtorus_addsub_params
prepare_addsub(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, UINT64 limb_precision) {
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
    mpn_add_n(dest.limbs + prep.doffset,
              a.limbs + prep.aoffset,
              b.limbs + prep.boffset,
              prep.limb_precision);
}

void sub_prep(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, const bigtorus_addsub_params &prep) {
    mpn_sub_n(dest.limbs + prep.doffset,
              a.limbs + prep.aoffset,
              b.limbs + prep.boffset,
              prep.limb_precision);
}

void add(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, UINT64 limb_precision) {
    add_prep(dest, a, b, prepare_addsub(dest, a, b, limb_precision));
}

void sub(BigTorusRef dest, const BigTorusRef &a, const BigTorusRef &b, UINT64 limb_precision) {
    sub_prep(dest, a, b, prepare_addsub(dest, a, b, limb_precision));
}

void zero(BigTorusRef dest) {
    mpn_zero(dest.limbs, dest.params->torus_limbs);
}

void subMul(BigTorusRef out, __int128 a, const BigTorusRef &in, UINT64 out_limb_prec) {
    if (out_limb_prec == NA) {
        out_limb_prec = out.params->torus_limbs;
    } else {
        assert_dramatically(out_limb_prec <= out.params->torus_limbs, "requested output precision too high");
    }
    assert_dramatically(out_limb_prec >= 1, "precision too low");
    const UINT64 in_limb_prec = out_limb_prec + 2; //input is 128-bit more precise
    assert_dramatically(in_limb_prec <= in.params->torus_limbs, "requested input precision too high");

    UINT64 *out_ptr = out.limbs + out.params->torus_limbs - out_limb_prec;
    UINT64 *in_ptr = in.limbs + in.params->torus_limbs - in_limb_prec;
    UINT64 *tmp = new UINT64[in_limb_prec];
    const UINT64 *a64 = (UINT64 *) &a; //a interpreted as two 64-bit unsigned ints
    if (a >= 0) {
        mpn_mul_1(tmp, in_ptr, in_limb_prec, a64[0]);
        mpn_addmul_1(tmp + 1, in_ptr, in_limb_prec - 1, a64[1]);
        mpn_sub_n(out_ptr, out_ptr, tmp + 2, out_limb_prec);
    } else {
        a = -a;
        mpn_mul_1(tmp, in_ptr, in_limb_prec, a64[0]);
        mpn_addmul_1(tmp + 1, in_ptr, in_limb_prec - 1, a64[1]);
        mpn_add_n(out_ptr, out_ptr, tmp + 2, out_limb_prec);
    }
    delete[] tmp;
}

void mul64(BigTorusRef out, int64_t a, const BigTorusRef &in, UINT64 limb_prec) {
    if (limb_prec == NA) {
        limb_prec = out.params->torus_limbs;
    } else {
        assert_dramatically(limb_prec <= out.params->torus_limbs, "requested output precision too high");
    }
    assert_dramatically(limb_prec >= 1, "precision too low");
    assert_dramatically(limb_prec <= in.params->torus_limbs, "requested input precision too high");

    UINT64 *out_ptr = out.limbs + out.params->torus_limbs - limb_prec;
    UINT64 *in_ptr = in.limbs + in.params->torus_limbs - limb_prec;
    if (a >= 0) {
        mpn_mul_1(out_ptr, in_ptr, limb_prec, a);
    } else {
        a = -a;
        mpn_mul_1(out_ptr, in_ptr, limb_prec, a);
        mpn_neg(out_ptr, out_ptr, limb_prec);
    }
}

