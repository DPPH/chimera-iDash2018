#include "TRLwe.h"

void pubKS(TRLwe &out, TLwe &in, pubKsKey &ks, const UINT64 out_prec_limbs) {
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
            subMul(out, aij, ks.kskey[i][j], out_prec_limbs);
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
    return 0;
}

void trivial(TRLwe &out, const BigTorusRef &in, const UINT64 out_limb_prec) {
    zero(out.a[0]);
    const_poly(out.a[1], in, out_limb_prec);
}

void subMul(TRLwe &out, __int128 aij, const TRLwe &in, const UINT64 out_limb_prec) {
    subMul(out.a[0], aij, in.a[0], out_limb_prec);
    subMul(out.a[1], aij, in.a[1], out_limb_prec);
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
}

void native_encrypt(TRLwe &reps, const BigTorusPolynomial &plaintext, const TLweKey &key, UINT64 alpha_bits) {

    const UINT64 N = reps.params.N;
    const UINT64 alpha_limbs = limb_precision(alpha_bits);

    BigTorusPolynomial &b = reps.a[1];
    // b = plaintext + sum s_i a_i
    copy(b, plaintext, alpha_limbs);
    random(reps.a[0], alpha_limbs);
    BigTorusPolynomial temp(N, alpha_limbs);
    int64_t *K = new int64_t[N];

    for (UINT64 i = 0; i < N; i++) {
        K[i] = key.key[i];
    }
    naive_external_product(temp, K, reps.a[0], alpha_limbs);
    add(b, b, temp, alpha_limbs);
    //randomize below bit alpha (noise)
    add_noise(b, alpha_bits, alpha_limbs);
    delete[] K;
}


void native_phase(BigTorusPolynomial &reps, const TRLwe &c, const TLweKey &key, UINT64 alpha_bits) {
    const UINT64 N = c.params.N;
    const UINT64 alpha_limbs = limb_precision(alpha_bits);

    const BigTorusPolynomial &b = c.a[1];
    copy(reps, b, alpha_limbs);

    BigTorusPolynomial temp(N, alpha_limbs);
    int64_t *K = new int64_t[N];


    for (UINT64 i = 0; i < N; i++) {
        K[i] = key.key[i];
    }

    naive_external_product(temp, K, c.a[0], alpha_limbs);
    sub(reps, b, temp, alpha_limbs);
    delete[] K;
}

void zero(TRLwe &out) {
    zero(out.a[0]);
    zero(out.a[1]);
}


std::shared_ptr<pubKsKey>
ks_keygen(const TRLweParams &out_params, const TLweParams &in_params, const TLweKey &in_key, const TLweKey &out_key,
          const UINT64 alpha_bits) {

    UINT64 ldec = (out_params.fixp_params.torus_limbs + 1) / 2;
    pubKsKey *reps = new pubKsKey(in_params, out_params, ldec);
    BigTorusPolynomial plaintext(out_params.N, out_params.fixp_params);

    for (UINT64 i = 0; i < in_params.N; i++) {
        for (UINT64 j = 1; j <= ldec; j++) {
            zero(plaintext);
            if (in_key.key[i] == 1) { plaintext.getAT(0).limbs_end[-2 * j] = 1; }
            native_encrypt(reps->kskey[i][j], plaintext, out_key, alpha_bits);
        }
    }

// sum Bg/2 Bg^-i
    for (UINT64 i = 1; i <= ldec; i++) {
        reps->bitDecomp_in_offset.limbs_end[-2 * i + 1] = 0x8000000000000000UL;
        reps->bitDecomp_in_offset.limbs_end[-2 * i] = 0;
    }

    reps->bitDecomp_out_offset = __int128(1) << 127; // -Bg/2
    return std::shared_ptr<pubKsKey>(reps);
}




pubKsKey::pubKsKey(const TLweParams &in_params, const TRLweParams &out_params, const UINT64 l_dec) :
        in_params(in_params),
        out_params(out_params),
        l_dec(l_dec),
        kskey(new_TRLweMatrix(in_params.N, l_dec, out_params)),
        bitDecomp_in_offset(in_params.fixp_params),
        bitDecomp_out_offset(0) {}

pubKsKey::~pubKsKey() {
    delete_TRLweMatrix(kskey, in_params.N, l_dec, out_params);

}
