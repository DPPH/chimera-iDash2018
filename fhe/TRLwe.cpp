#include <cassert>
#include <vector>
#include "TRLwe.h"


// see if this out_prec_limbs still make sense?
void pubKS128(TRLwe &out, const TLwe &in, const pubKsKey128 &ks, const UINT64 out_prec_limbs) {
    //const TRLweParams& out_params = out.params;
    const TLweParams &in_params = in.params;
    //const int64_t ks_prec_limbs = out_prec_limbs+2; //ks has 128-bit more precision
    const BigTorus &bitDecomp_in_offset = ks.bitDecomp_in_offset; // sum Bg/2 Bg^i
    __int128 bitDecomp_out_offset = ks.bitDecomp_out_offset; // -Bg/2

    BigTorus tmpDec(in_params.fixp_params); //temp variable

    // out = trivial(b)
    trivial(out, in.getBT(), out_prec_limbs);
    for (UINT64 i = 0; i < in_params.N; i++) {
        add(tmpDec, in.getAT(i), bitDecomp_in_offset);
        for (UINT64 j = 1; j <= ks.l_dec; ++j) {
            // coef aij of the decomposition
            __int128 aij = bitdecomp_coef128(tmpDec, j, out_prec_limbs) + bitDecomp_out_offset;
            // out = out - aij . ks_ij
            subMul128(out, aij, ks.kskey[i][j - 1], out_prec_limbs);
        }
    }
}

// see if this out_prec_limbs still make sense?
void pubKS32(TRLwe &out, const TLwe &in, const pubKsKey32 &ks, const UINT64 out_prec_limbs) {
    //const TRLweParams& out_params = out.params;
    const TLweParams &in_params = in.params;
    //const int64_t ks_prec_limbs = out_prec_limbs+2; //ks has 128-bit more precision
    //const BigTorus &bitDecomp_in_offset = ks.bitDecomp_in_offset; // sum Bg/2 Bg^i
    //int64_t bitDecomp_out_offset = ks.bitDecomp_out_offset; // -Bg/2

    BigTorus tmpDec(in_params.fixp_params); //temp variable

    // out = trivial(b)
    trivial(out, in.getBT(), out_prec_limbs);
    for (UINT64 i = 0; i < in_params.N; i++) {
        bitdecomp_signed_offset32_apply(tmpDec, in.getAT(i));
        for (UINT64 j = 1; j <= ks.l_dec; ++j) {
            // coef aij of the decomposition
            //int64_t aij = bitdecomp_coef32(tmpDec, j, out_prec_limbs) + bitDecomp_out_offset;
            int64_t aij = bitdecomp_signed_coef32(tmpDec, j, out_prec_limbs);
            // out = out - aij . ks_ij
            subMul64(out, aij, ks.kskey[i][j - 1], out_prec_limbs);
        }
    }
}


__int128 bitdecomp_coef128(const BigTorusRef &tmpDec, UINT64 j, const UINT64 limb_prec) {
    const UINT64 nlimbs = tmpDec.params.torus_limbs;
    if (2 * j <= nlimbs) {
        return *(__int128 *) (tmpDec.limbs_end - 2 * j);
    } else if (2 * j == nlimbs + 1) {
        return ((__int128) *(tmpDec.limbs_end - nlimbs)) << (__int128(64));
    }
    //fill with zero beyond limb_prec
    return 0;
}

void bitdecomp_signed_offset32_apply(BigTorusRef reps, const BigTorusRef &source) {
    int64_t dlimb = reps.params.torus_limbs;
    int64_t slimb = source.params.torus_limbs;
    static const UINT64 LIMB_OFFSET = 0x8000000080000000UL; //2^31 + 2^63
    static std::vector<UINT64> offsetV;
    if (offsetV.size() < UINT64(dlimb)) {
        offsetV.resize(dlimb, LIMB_OFFSET);
    }
    //copy source to destination, filling with zeros
    if (dlimb <= slimb) {
        for (int64_t i = 1; i <= dlimb; i++) {
            reps.limbs_end[-i] = source.limbs_end[-i];
        }
    } else {
        for (int64_t i = 1; i <= slimb; i++) {
            reps.limbs_end[-i] = source.limbs_end[-i];
        }
        for (int64_t i = slimb + 1; i <= dlimb; i++) {
            reps.limbs_end[-i] = 0;
        }
    }
    //add offset
    mpn_add_n(reps.limbs_end - dlimb, reps.limbs_end - dlimb, offsetV.data(), dlimb);
}


int64_t bitdecomp_coef32(const BigTorusRef &tmpDec, UINT64 j, const UINT64 limb_prec) {
    assert(j >= 1);
    const UINT64 nlimbs = tmpDec.params.torus_limbs;
    if (j <= 2 * nlimbs) {
        if (j % 2 == 0) {
            return int64_t(tmpDec.limbs_end[-j / 2] & 0xFFFFFFFFUL);
        } else {
            return int64_t(tmpDec.limbs_end[-(j + 1) / 2] >> 32UL);
        }
    } else return 0;
}

int64_t bitdecomp_signed_coef32(const BigTorusRef &tmpDec, UINT64 j, const UINT64 limb_prec) {
    static const int64_t bitDecomp_out_offset = -(int64_t(1) << int64_t(31)); // -Bg/2
    return bitdecomp_coef32(tmpDec, j, limb_prec) + bitDecomp_out_offset;
}

void trivial(TRLwe &out, const BigTorusRef &in, const UINT64 out_limb_prec) {
    zero(out.a[0]);
    const_poly(out.a[1], in, out_limb_prec);
}

void subMul128(TRLwe &out, __int128 aij, const TRLwe &in, const UINT64 out_limb_prec) {
    subMul128(out.a[0], aij, in.a[0], out_limb_prec);
    subMul128(out.a[1], aij, in.a[1], out_limb_prec);
}

void subMul64(TRLwe &out, int64_t aij, const TRLwe &in, const UINT64 out_limb_prec) {
    subMul64(out.a[0], aij, in.a[0], out_limb_prec);
    subMul64(out.a[1], aij, in.a[1], out_limb_prec);
}

TRLweParams::TRLweParams(const UINT64 N, const BigTorusParams &fixp_params) :
        TLweParams(N, fixp_params) {}

TRLwe::TRLwe(const TRLweParams &params) :
        params(params),
        a((BigTorusPolynomial *) malloc(2 * sizeof(BigTorusPolynomial))) {
    new(a) BigTorusPolynomial(params.N, params.fixp_params);
    new(a + 1) BigTorusPolynomial(params.N, params.fixp_params);
}

TRLwe::~TRLwe() {
    a->~BigTorusPolynomial();
    (a + 1)->~BigTorusPolynomial();
    free(a);
}

void native_encrypt(TRLwe &reps, const BigTorusPolynomial &plaintext, const TLweKey &key, UINT64 alpha_bits) {

    const UINT64 N = reps.params.N;
    const UINT64 alpha_limbs = limb_precision(alpha_bits);
    const BigTorusParams bt_out_params(alpha_limbs);

    BigTorusPolynomial &b = reps.a[1];
    // b = plaintext + sum s_i a_i
    copy(b, plaintext, alpha_limbs);
    random(reps.a[0], alpha_limbs);
    BigTorusPolynomial temp(N, bt_out_params);
    int64_t *K = new int64_t[N];

    for (UINT64 i = 0; i < N; i++) {
        K[i] = key.key[i];
    }
    fft_external_product(temp, K, reps.a[0], 1, alpha_limbs);
    add(b, b, temp, alpha_limbs);
    //randomize below bit alpha (noise)
    add_noise(b, alpha_bits, alpha_limbs);

    delete[] K;
}


void native_phase(BigTorusPolynomial &reps, const TRLwe &c, const TLweKey &key, UINT64 alpha_bits) {
    const UINT64 N = c.params.N;
    const UINT64 alpha_limbs = limb_precision(alpha_bits);
    BigTorusParams bt_params(alpha_limbs);

    const BigTorusPolynomial &b = c.a[1];
    copy(reps, b, alpha_limbs);

    BigTorusPolynomial temp(N, bt_params);
    int64_t *K = new int64_t[N];


    for (UINT64 i = 0; i < N; i++) {
        K[i] = key.key[i];
    }

    fft_external_product(temp, K, c.a[0], 1, alpha_limbs);
    sub(reps, b, temp, alpha_limbs);
    delete[] K;
}

void zero(TRLwe &out) {
    zero(out.a[0]);
    zero(out.a[1]);
}

TRLwe **new_TRLweMatrix(UINT64 rows, UINT64 cols, const TRLweParams &params) {
    uint8_t *r = (uint8_t *) malloc(rows * sizeof(TRLwe *) + rows * cols * sizeof(TRLwe));
    TRLwe **reps = (TRLwe **) r;
    TRLwe *reps_raw = (TRLwe *) (r + rows * sizeof(TRLwe *));
    for (UINT64 i = 0; i < rows * cols; i++) {
        new(reps_raw + i) TRLwe(params);
    }
    for (UINT64 j = 0; j < rows; j++) {
        reps[j] = reps_raw + cols * j;
    }
    return reps;
}

void delete_TRLweMatrix(TRLwe **data, UINT64 rows, UINT64 cols, const TRLweParams &params) {
    uint8_t *r = (uint8_t *) data;
    TRLwe *reps_raw = (TRLwe *) (r + rows * sizeof(TRLwe *));
    for (UINT64 i = 0; i < rows * cols; i++) {
        (reps_raw + i)->~TRLwe();
    }
    free(r);
}


std::shared_ptr<pubKsKey128>
ks_keygen128(const TRLweParams &out_params, const TLweParams &in_params,
             const TLweKey &in_key, const TLweKey &out_key,
             const UINT64 out_alpha_bits) {
    //automatic parameter deduction
    const double out_variance_bits = out_alpha_bits * 2;
    const double in_variance_bits = out_variance_bits + 1;     // 1/2 from the input noise
    const double errdec_variance_bits = out_variance_bits + 2; // 1/4 from decomp error
    const double kssum_variance_bits = out_variance_bits + 2;  // 1/4 from the big sum


    //smallest decomposition in base 2^128 giving an error variance <= errdec_variance_bits
    const UINT64 ldec = ceil(errdec_variance_bits / 2. / 128.);

    //the input noise variance must be smaller than in_variance_bits
    // -> input nblimbs must be larger than limb_prec(variance/2)
    assert_dramatically(in_params.fixp_params.torus_limbs >= UINT64(limb_precision(in_variance_bits / 2)));

    //smallest variance for the keyswitch key
    const double ks_variance_bits = kssum_variance_bits + log2(in_params.N) + log2(ldec) + 2 * 127.;
    const UINT64 ks_alpha_bits = UINT64(ks_variance_bits / 2);
    const UINT64 ks_nblimbs = UINT64(limb_precision(ks_alpha_bits) + 1) & UINT64(-2); //must be even
    TRLweParams ks_params(out_params.N, ks_nblimbs);

    //the output nblimbs must be larger than limb_prec(out_alpha_bits)
    assert_dramatically(out_params.fixp_params.torus_limbs >= UINT64(limb_precision(out_alpha_bits)));

    pubKsKey128 *reps = new pubKsKey128(in_params, out_params, TRLweParams(ks_params), ldec);

    //plaintext must have the same precision as the keyswitch key
    BigTorusPolynomial plaintext(ks_params.N, ks_params.fixp_params);
    for (UINT64 i = 0; i < in_params.N; i++) {
        for (UINT64 j = 1; j <= ldec; j++) {
            // plaintext = ks[i] / Bg^j
            zero(plaintext);
            if (in_key.key[i] == 1) { plaintext.getAT(0).limbs_end[-2 * j] = 1; }
            native_encrypt(reps->kskey[i][j - 1], plaintext, out_key, ks_alpha_bits);
        }
    }

    // sum Bg/2 Bg^-i
    zero(reps->bitDecomp_in_offset);
    for (UINT64 i = 1; i <= reps->in_params.fixp_params.torus_limbs; i += 2) {
        reps->bitDecomp_in_offset.limbs_end[-i] = 0x8000000000000000UL;
    }

    reps->bitDecomp_out_offset = (__int128(1) << __int128(127)); // -Bg/2
    return std::shared_ptr<pubKsKey128>(reps);
}

std::shared_ptr<pubKsKey32>
ks_keygen32(const TRLweParams &out_params, const TLweParams &in_params, const TLweKey &in_key, const TLweKey &out_key,
            const UINT64 out_alpha_bits) {
    //automatic parameter deduction
    const double out_variance_bits = out_alpha_bits * 2;
    const double in_variance_bits = out_variance_bits + 1;     // 1/2 from the input noise
    const double errdec_variance_bits = out_variance_bits + 2; // 1/4 from decomp error
    const double kssum_variance_bits = out_variance_bits + 2;  // 1/4 from the big sum


    //smallest decomposition in base 2^32 giving an error variance <= errdec_variance_bits
    const UINT64 ldec = ceil(errdec_variance_bits / 2. / 32.);

    //the input noise variance must be smaller than in_variance_bits
    // -> input nblimbs must be larger than limb_prec(variance/2)
    assert_dramatically(in_params.fixp_params.torus_limbs >= UINT64(limb_precision(in_variance_bits / 2)));

    //smallest variance for the keyswitch key
    const double ks_variance_bits = kssum_variance_bits + log2(in_params.N) + log2(ldec) + 2 * 31.;
    const UINT64 ks_alpha_bits = UINT64(ks_variance_bits / 2);
    const UINT64 ks_nblimbs = limb_precision(ks_alpha_bits);
    TRLweParams ks_params(out_params.N, ks_nblimbs);

    //the output nblimbs must be larger than limb_prec(out_alpha_bits)
    assert_dramatically(out_params.fixp_params.torus_limbs >= UINT64(limb_precision(out_alpha_bits)));

    pubKsKey32 *reps = new pubKsKey32(in_params, out_params, TRLweParams(ks_params), ldec);

    //plaintext must have the same precision as the keyswitch key
    BigTorusPolynomial plaintext(ks_params.N, ks_params.fixp_params);
    for (UINT64 i = 0; i < in_params.N; i++) {
        for (UINT64 j = 1; j <= ldec; j++) {
            // plaintext = ks[i] / Bg^j
            zero(plaintext);
            if (in_key.key[i] == 1) {
                if (j % 2 == 0) {
                    plaintext.getAT(0).limbs_end[-j / 2] = 1;
                } else {
                    plaintext.getAT(0).limbs_end[-(j + 1) / 2] = (1ul << 32ul);
                }
            }
            native_encrypt(reps->kskey[i][j - 1], plaintext, out_key, ks_alpha_bits);
        }
    }

    // sum Bg/2 Bg^-i
    zero(reps->bitDecomp_in_offset);
    for (UINT64 i = 1; i <= reps->in_params.fixp_params.torus_limbs; i++) {
        reps->bitDecomp_in_offset.limbs_end[-i] = 0x8000000080000000UL;
    }

    reps->bitDecomp_out_offset = -(int64_t(1) << int64_t(31)); // -Bg/2
    return std::shared_ptr<pubKsKey32>(reps);
}


pubKsKey128::pubKsKey128(
        const TLweParams &in_params,
        const TRLweParams &out_params,
        TRLweParams &&ks_params,
        const UINT64 l_dec) :
        in_params(in_params),
        out_params(out_params),
        ks_params(ks_params),
        l_dec(l_dec),
        kskey(new_TRLweMatrix(in_params.N, l_dec, this->ks_params)), //warning: this-> is essential here!
        bitDecomp_in_offset(in_params.fixp_params),
        bitDecomp_out_offset(0) {}

pubKsKey128::~pubKsKey128() {
    delete_TRLweMatrix(kskey, in_params.N, l_dec, out_params);
}

pubKsKey32::pubKsKey32(const TLweParams &in_params, const TRLweParams &out_params, TRLweParams &&ks_params,
                       const UINT64 l_dec) :
        in_params(in_params),
        out_params(out_params),
        ks_params(ks_params),
        l_dec(l_dec),
        kskey(new_TRLweMatrix(in_params.N, l_dec, this->ks_params)), //warning: this-> is essential here!
        bitDecomp_in_offset(in_params.fixp_params),
        bitDecomp_out_offset(0) {}

pubKsKey32::~pubKsKey32() {
    delete_TRLweMatrix(kskey, in_params.N, l_dec, out_params);
}
