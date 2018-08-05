#include "TLwe.h"
#include "arithmetic.h"

TLweParams::TLweParams(const uint64_t N, const BigFixPParams &fixp_params) : N(N), fixp_params(fixp_params) {}

TLwe::TLwe(const TLweParams &params) :
        limbs(new uint64_t[(params.N + 1) * params.fixp_params.torus_limbs * sizeof(uint64_t)]),
        params(params) {
}

TLwe::~TLwe() {
    delete[] limbs;
}

BigTorusRef TLwe::getAT(uint64_t i) {
    return BigTorusRef(limbs + i * params.fixp_params.torus_limbs, &params.fixp_params);
}

BigTorusRef TLwe::getAT(uint64_t i) const {
    return BigTorusRef(limbs + i * params.fixp_params.torus_limbs, &params.fixp_params);
}

BigTorusRef TLwe::getBT() {
    return BigTorusRef(limbs + params.N * params.fixp_params.torus_limbs,
                       &params.fixp_params);
}

BigTorusRef TLwe::getBT() const {
    return BigTorusRef(limbs + params.N * params.fixp_params.torus_limbs,
                       &params.fixp_params);
}

TLweKey::TLweKey(const TLweParams &params) :
        key(new int8_t[params.N]),
        params(params) {
}

TLweKey::~TLweKey() {
    delete[] key;
}

std::shared_ptr<TLweKey> tlwe_keygen(const TLweParams &params) {
    TLweKey *reps = new TLweKey(params);
    const uint64_t N = params.N;
    for (uint64_t i = 0; i < N; i++) {
        reps->key[i] = random_bit();
    }
    return std::shared_ptr<TLweKey>(reps);
}

void zero(TLwe &tlwe) {
    const uint64_t Np = tlwe.params.N + 1;
    const uint64_t limbSize = tlwe.params.fixp_params.torus_limbs;
    mpn_zero(tlwe.limbs, Np * limbSize);
}

void native_encrypt(TLwe &reps, const BigTorusRef &plaintext, const TLweKey &key, uint64_t alpha_bits) {
    const uint64_t N = reps.params.N;
    const uint64_t alpha_limbs = limb_precision(alpha_bits);
    auto b = reps.getBT();
    // b = plaintext + sum s_i a_i
    copy(b, plaintext, alpha_limbs);
    auto prep = prepare_addsub(b, b, b, alpha_limbs);
    for (uint64_t i = 0; i < N; i++) {
        auto ai = reps.getAT(i);
        random(ai, alpha_limbs);
        if (key.key[i]) {
            add_prep(b, b, ai, prep);
        }
    }
    //randomize below bit alpha (noise)
    add_noise(b, alpha_bits, alpha_limbs);
}

void native_phase(BigTorusRef reps, const TLwe &tlwe, const TLweKey &key, uint64_t alpha_bits) {
    const uint64_t N = key.params.N;
    const uint64_t alpha_limbs = limb_precision(alpha_bits);
    auto b = tlwe.getBT();
    copy(reps, b, alpha_limbs);
    auto prep = prepare_addsub(reps, reps, b, alpha_limbs);
    for (uint64_t i = 0; i < N; i++) {
        auto ai = tlwe.getAT(i);
        if (key.key[i]) {
            sub_prep(reps, reps, ai, prep);
        }
    }
}

void slot_encrypt(TLwe &reps, const NTL::RR &plaintext, const TLweKey &key, uint64_t alpha_bits) {
    BigFixP tmp(&reps.params.fixp_params);
    to_fixP(tmp, plaintext);
    native_encrypt(reps, BigTorusRef(tmp.limbs_raw, tmp.params), key, alpha_bits);
}

NTL::RR slot_decrypt(const TLwe &tlwe, const TLweKey &key, uint64_t alpha_bits) {
    BigFixP tmp(&tlwe.params.fixp_params);
    native_phase(BigTorusRef(tmp.limbs_raw, tmp.params), tlwe, key, alpha_bits);
    return to_RR(tmp);
}
