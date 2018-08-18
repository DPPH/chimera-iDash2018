#include <cassert>
#include "BigReal.h"
#include "BigFFT.h"

BigReal::BigReal(UINT64 nblimbs) : nblimbs(nblimbs) {
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

void div2ui(BigReal &dest, const BigReal &a, UINT64 b) {
    mpz_cdiv_q_2exp(dest.value, a.value, b);
}

void to_BigReal(BigReal &dest, const BigTorusRef &v) {
    const int dsize = dest.nblimbs;
    const int vsize = v.params.torus_limbs;
    const int vcopysize = std::min(dsize, vsize);
    const int zerofill = dsize - vcopysize;
    //copy the nblimbs most significant limbs
    mpz_realloc2(dest.value, dsize * BITS_PER_LIMBS);
    assert(dest.value->_mp_alloc >= dsize);
    //dest.value->_mp_d = (mp_limb_t *) malloc(dsize * BYTES_PER_LIMBS);
    //fill the zeros if dest is larger
    if (zerofill > 0) {
        mpn_zero(dest.value->_mp_d, zerofill);
    }
    //negate or copy limbs from v
    bool vnegative = v.limbs_end[-1] >> 63u;
    if (vnegative) {
        //v is negative
        mpn_neg(dest.value->_mp_d + zerofill, v.limbs_end - vcopysize, vcopysize);
    } else {
        //v is positive
        mpn_copyi(dest.value->_mp_d + zerofill, v.limbs_end - vcopysize, vcopysize);
    }
    //find the right size (for the mpz_t)
    int dl;
    for (dl = dsize; dl >= 1; dl--) {
        if (dest.value->_mp_d[dl - 1] != 0) break;
    }
    dest.value->_mp_size = vnegative ? -dl : dl;
}

void to_BigReal(BigReal &dest, int64_t a, UINT64 a_nbits) {
    assert(a_nbits <= 64 && a_nbits >= 1);
    if (a == 0) {
        mpz_set_ui(dest.value, 0);
        return;
    }
    const int dsize = dest.nblimbs;
    mpz_realloc2(dest.value, dsize * BITS_PER_LIMBS);
    assert(dest.value->_mp_alloc >= dsize);
    mpn_zero(dest.value->_mp_d, dsize);
    if (a >= 0) {
        dest.value->_mp_d[dsize - 1] = (UINT64(a) << UINT64(64 - a_nbits));
        dest.value->_mp_size = dsize;
    } else {
        a = -a;
        dest.value->_mp_d[dsize - 1] = (UINT64(a) << UINT64(64 - a_nbits));
        dest.value->_mp_size = -dsize;
    }
}


void copy(BigReal &dest, const BigReal &a) {
    mpz_set(dest.value, a.value);
}

BigReal *new_BigReal_array(UINT64 n, UINT64 nblimbs) {
    BigReal *reps = (BigReal *) malloc(n * sizeof(BigReal));
    for (UINT64 i = 0; i < n; i++) {
        new(reps + i) BigReal(nblimbs);
    }
    return reps;
}

void delete_BigReal_array(UINT64 n, BigReal *array) {
    for (UINT64 i = 0; i < n; i++) {
        (array + i)->~BigReal();
    }
    free(array);
}

#include <NTL/ZZ_limbs.h>
#include <gmp.h>
#include <vector>

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
    return to_RR(vv) / pow(to_RR(2), to_RR(long(v.nblimbs * BITS_PER_LIMBS)));
}

void to_BigReal(BigReal &dest, const NTL::RR &v) {
    assert_dramatically(abs(v) <= 1, "Bad range for BigReal");
    RR::SetPrecision(dest.nblimbs * BITS_PER_LIMBS);
    ZZ vv = RoundToZZ(v * pow(to_RR(2), to_RR(long(dest.nblimbs * BITS_PER_LIMBS))));
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

void to_BigTorus(BigTorusRef dest, const BigReal &a, UINT64 lshift_bits, UINT64 out_precision_limbs) {
    assert(lshift_bits >= 0 && lshift_bits < 64);
    if (out_precision_limbs == NA) {
        out_precision_limbs = dest.params.torus_limbs;
    } else {
        assert_dramatically(dest.params.torus_limbs >= out_precision_limbs);
    }
    UINT64 *tmp = new UINT64[out_precision_limbs + 1];
    //copy all limbs from a to tmp
    //a.value->_mp_d[a.nblimbs-1] -> tmp[out_precision_limbs]
    //a.value->_mp_d[a.nblimbs-2] -> tmp[out_precision_limbs-1]
    //...
    int64_t asize = abs(a.value->_mp_size);
    for (UINT64 i = 0; i <= out_precision_limbs; i++) {
        if (int64_t(a.nblimbs - i - 1) < asize && int64_t(a.nblimbs - i - 1) >= 0)
            tmp[out_precision_limbs - i] = a.value->_mp_d[a.nblimbs - 1 - i];
        else
            tmp[out_precision_limbs - i] = 0;
    }
    //negate if a is negative
    if (a.value->_mp_size < 0) {
        mpn_neg(tmp, tmp, out_precision_limbs + 1);
    }
    //do the shift
    if (lshift_bits > 0)
        mpn_lshift(tmp, tmp, out_precision_limbs + 1, lshift_bits);
    //copy the final limbs
    mpn_copyi(dest.limbs_end - out_precision_limbs, tmp + 1, out_precision_limbs);
    //free resources
    delete[] tmp;
}

void zero(BigReal &dest) {
    mpz_set_ui(dest.value, 0);
}


void precise_conv_toBigReal(BigReal &reps, const BigTorusRef &b, int64_t lshift, int64_t msbToKeep) {
    int64_t blimbs = b.params.torus_limbs;
    int64_t shift_mod_64 = lshift % BITS_PER_LIMBS;
    int64_t shift_over_64 = lshift / BITS_PER_LIMBS;
    int64_t newblimbs = blimbs - shift_over_64; //new size of b
    if (newblimbs <= 0) {
        //corner case: the answer is 0
        zero(reps);
        return;
    }
    //shift b by lshift positions and store the limbs in btmp
    std::vector<UINT64> btmp(newblimbs);
    UINT64 *newbbegin = btmp.data();
    UINT64 *newbend = newbbegin + newblimbs;
    mpn_copyi(newbbegin, b.limbs_end - blimbs, newblimbs);
    if (shift_mod_64 != 0) {
        mpn_lshift(newbbegin, newbbegin, newblimbs, shift_mod_64);
    }

    int64_t nlimbs_to_copy = limb_precision(msbToKeep);
    assert_dramatically(int64_t(reps.nblimbs) >= nlimbs_to_copy, "a does is not large enough to store the result");
    mpz_realloc2(reps.value, reps.nblimbs * BITS_PER_LIMBS);
    assert(int64_t(reps.value->_mp_alloc) >= int64_t(reps.nblimbs));
    UINT64 *dest_end = reps.value->_mp_d + reps.nblimbs;
    //copy the msb limbs
    for (int64_t i = 1; i <= nlimbs_to_copy; i++) {
        if (i <= int64_t(newblimbs))
            dest_end[-i] = newbend[-i];
        else
            dest_end[-i] = 0;

    }
    //fill the lsb with zeros if needed
    for (int64_t i = nlimbs_to_copy + 1; i <= int64_t(reps.nblimbs); i++) {
        dest_end[-i] = 0;
    }
    if (msbToKeep % 64 != 0) {
        //zeroify 64-(msb%64) bits on the least significant limb copied
        int64_t k = 64 - (msbToKeep % 64);
        int64_t mask = -(1l << k);
        //add 1<<(k-1)
        mpn_add_1(dest_end - nlimbs_to_copy, dest_end - nlimbs_to_copy, nlimbs_to_copy, (1ul << uint64_t(k - 1)));
        //then zeroify
        dest_end[-nlimbs_to_copy] &= mask;
    } else {
        //recenter too
        if (nlimbs_to_copy < newblimbs && (newbend[-nlimbs_to_copy - 1] & 0x8000000000000000ul) != 0) {
            mpn_add_1(dest_end - nlimbs_to_copy, dest_end - nlimbs_to_copy, nlimbs_to_copy, 1ul);
        }
    }
    //if the result is negative, negate it
    bool vneg = ((dest_end[-1] & 0x8000000000000000UL) != 0);
    if (vneg) {
        mpn_neg(dest_end - nlimbs_to_copy, dest_end - nlimbs_to_copy, nlimbs_to_copy);
    }
    //re-compute the actual size of a
    int64_t i;
    for (i = reps.nblimbs; i >= 1; i--) {
        if (reps.value->_mp_d[i - 1] != 0) break;
    }
    reps.value->_mp_size = vneg ? -i : i;
}

void fft_BigRealPoly_product(BigReal *reps, BigReal *a, BigReal *b, int64_t N, int64_t fft_nlimbs) {

    int64_t n = 2 * N;
    int64_t Ns2 = N / 2;

    BigComplex *ca = new_BigComplex_array(Ns2, fft_nlimbs);
    BigComplex *cb = new_BigComplex_array(Ns2, fft_nlimbs);

    const BigComplex *powomega = fftAutoPrecomp.omega(n, fft_nlimbs);
    const BigComplex *powombar = fftAutoPrecomp.omegabar(n, fft_nlimbs);
    iFFT(ca, a, n, powomega);
    iFFT(cb, b, n, powomega);

    for (int i = 0; i < Ns2; i++) {
        mulTo(ca[i], cb[i]);
    }

    FFT(reps, ca, n, powombar);
    delete_BigComplex_array(Ns2, ca);
    delete_BigComplex_array(Ns2, cb);
}

void BigRealPoly_addTo(BigReal *reps, BigReal *in, int64_t N, int64_t fft_nlimbs) {

    for (int i = 0; i < N; i++) {
        add(reps[i], in[i], reps[i]);
    }
}

void shift_toBigTorus(BigTorusRef out, const BigReal &a, int64_t left_shift) {
    //treat the sero corner case first
    if (a.value->_mp_size == 0) {
        zero(out);
        return;
    }
    // out = a * 2^left_shift / 2^64.(a.nblimbs)
    const int64_t outlimbs = out.params.torus_limbs;
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mul_2exp(tmp, a.value, left_shift);
    const int64_t tlimbs = abs(tmp->_mp_size);
    //copy  outlimbs limbs starting from position alimbs
    for (int64_t i = 1; i <= outlimbs; i++) {
        int64_t j = a.nblimbs - i;
        if (j >= tlimbs) out.limbs_end[-i] = 0;
        else if (j >= 0) out.limbs_end[-i] = tmp->_mp_d[j];
        else out.limbs_end[-i] = 0;
    }
    //negate if a is negative
    if (a.value->_mp_size < 0) {
        mpn_neg(out.limbs_end - outlimbs, out.limbs_end - outlimbs, outlimbs);
    }
    //free resources
    mpz_clear(tmp);
}

