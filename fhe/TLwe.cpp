#include <gmp.h>
#include <cassert>
#include "TLwe.h"
#include "arithmetic.h"

using namespace std;

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

TLweKey::TLweKey(const UINT64 N) : key(new int64_t[N]) {

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

void serializeTLweParams(std::ostream &out, const TLweParams &tLweParams) {
    int64_t x = TLWE_PARAMS_SERIAL_ID;
    ostream_write_binary(out, &x, sizeof(int64_t));
    ostream_write_binary(out, &tLweParams.N, sizeof(UINT64));
    serializeBigTorusParams(out, tLweParams.fixp_params);
}

std::shared_ptr<TLweParams> deserializeTLweParams(std::istream &in) {
    int64_t magic;
    istream_read_binary(in, &magic, sizeof(int64_t));
    assert_dramatically(magic == TLWE_PARAMS_SERIAL_ID);
    UINT64 N;
    istream_read_binary(in, &N, sizeof(UINT64));
    shared_ptr<BigTorusParams> params = deserializeBigTorusParams(in);
    store_forever(params);
    TLweParams *tLweParams = new TLweParams(N, *params);
    return std::shared_ptr<TLweParams>(tLweParams);
}

void serializeTLweKey(std::ostream &out, const TLweKey &key, const UINT64 N) {
    int64_t x = TLWE_KEY_SERIAL_ID;
    ostream_write_binary(out, &x, sizeof(int64_t));
    for (int64_t i = 0; i < int64_t(N); i++) {
        int8_t coef = key.key[i];
        ostream_write_binary(out, &coef, sizeof(int8_t));
    }
}

std::shared_ptr<TLweKey> deserializeTLweKey(std::istream &in, const UINT64 N) {

    int64_t magic;
    istream_read_binary(in, &magic, sizeof(int64_t));
    assert_dramatically(magic == TLWE_KEY_SERIAL_ID);

    TLweKey *key = new TLweKey(N);

    for (int64_t i = 0; i < int64_t(N); i++) {
        int8_t coef;
        istream_read_binary(in, &coef, sizeof(int8_t));
        key->key[i] = coef;
    }
    return std::shared_ptr<TLweKey>(key);
}

void serializeTLwe(std::ostream &out, const TLwe &value) {
    int64_t x = TLWE_SERIAL_ID;
    ostream_write_binary(out, &x, sizeof(int64_t));
    serializeTLweParams(out, value.params);
    serializeBigTorusVectorContent(out, value);
}

std::shared_ptr<TLwe> deserializeTLwe(std::istream &in) {

    int64_t magic;
    istream_read_binary(in, &magic, sizeof(int64_t));
    assert_dramatically(magic == TLWE_SERIAL_ID);
    shared_ptr<TLweParams> params = deserializeTLweParams(in);
    store_forever(params);
    TLwe *reps = new TLwe(*params);
    deserializeBigTorusVectorContent(in, *reps);
    return shared_ptr<TLwe>(reps);
}
