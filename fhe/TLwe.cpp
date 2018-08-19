#include <gmp.h>
#include <cassert>
#include "TLwe.h"
#include "arithmetic.h"

TLweParams::TLweParams(const UINT64 N, const BigTorusParams &fixp_params) : N(N), fixp_params(fixp_params) {}

TLwe::TLwe(const TLweParams &params) :
        BigTorusVector(params.N + 1, params.fixp_params),
        params(params) {
}

TLwe::~TLwe() {
}

BigTorusRef TLwe::getBT() {
    return getAT(params.N);
}

BigTorusRef TLwe::getBT() const {
    return getAT(params.N);
}

TLweKey::TLweKey(const TLweParams &params) :
        key(new int64_t[params.N]) {
}

TLweKey::~TLweKey() {
    delete[] key;
}

std::shared_ptr<TLweKey> tlwe_keygen(const TLweParams &params) {
    TLweKey *reps = new TLweKey(params);
    const UINT64 N = params.N;
    for (UINT64 i = 0; i < N; i++) {
        reps->key[i] = random_bit();
    }
    return std::shared_ptr<TLweKey>(reps);
}

void zero(TLwe &tlwe) {
    const UINT64 Np = tlwe.params.N + 1;
    const UINT64 limbSize = tlwe.params.fixp_params.torus_limbs;
    mpn_zero(tlwe.limbs_end - limbSize, Np * limbSize);
}

void native_encrypt(TLwe &reps, const BigTorusRef &plaintext, const TLweKey &key, UINT64 alpha_bits) {
    const UINT64 N = reps.params.N;
    const UINT64 alpha_limbs = limb_precision(alpha_bits);
    auto b = reps.getBT();
    // b = plaintext + sum s_i a_i
    copy(b, plaintext, alpha_limbs);
    auto prep = prepare_addsub(b, b, b, alpha_limbs);
    for (UINT64 i = 0; i < N; i++) {
        auto ai = reps.getAT(i);
        random(ai, alpha_limbs);
        if (key.key[i]) {
            add_prep(b, b, ai, prep);
        }
    }
    //randomize below bit alpha (noise)
    add_noise(b, alpha_bits, alpha_limbs);
}

void native_phase(BigTorusRef reps, const TLwe &tlwe, const TLweKey &key, UINT64 alpha_bits) {
    const UINT64 N = tlwe.params.N;
    const UINT64 alpha_limbs = limb_precision(alpha_bits);
    auto b = tlwe.getBT();
    copy(reps, b, alpha_limbs);
    auto prep = prepare_addsub(reps, reps, b, alpha_limbs);
    for (UINT64 i = 0; i < N; i++) {
        auto ai = tlwe.getAT(i);
        if (key.key[i]) {
            sub_prep(reps, reps, ai, prep);
        }
    }
}

void slot_encrypt(TLwe &reps, const NTL::RR &plaintext, const TLweKey &key, UINT64 plaintext_precision) {
    assert(reps.params.fixp_params.level_expo > 0); //"level exponent is not set";
    int64_t alpha_bits = reps.params.fixp_params.level_expo + plaintext_precision;
    assert(int64_t(reps.params.fixp_params.torus_limbs) >=
           limb_precision(alpha_bits)); //"the bigtorus is not precise enough";
    BigTorus tmp(reps.params.fixp_params);
    to_fixP(tmp, plaintext);
    native_encrypt(reps, tmp, key, alpha_bits);
}

NTL::RR slot_decrypt(const TLwe &tlwe, const TLweKey &key) {
    assert(tlwe.params.fixp_params.level_expo > 0); //"level exponent is not set";
    BigTorus tmp(tlwe.params.fixp_params);
    native_phase(tmp, tlwe, key);
    return fixp_to_RR(tmp);
}
