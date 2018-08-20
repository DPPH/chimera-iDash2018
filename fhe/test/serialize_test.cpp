#include <include/gtest/gtest.h>
#include "../BigTorus.h"
#include "../BigTorusVector.h"
#include "../TLwe.h"
#include "../TRLwe.h"
#include <sstream>

using namespace std;

TEST(SERIALIZE_TEST, BigTorusParams) {
    BigTorusParams params(random(), random(), random());
    //serialize it into a string (instead of a file)
    ostringstream oss;
    serializeBigTorusParams(oss, params);
    string serial_result = oss.str();
    //deserialize
    istringstream iss(serial_result);
    shared_ptr<BigTorusParams> res_params = deserializeBigTorusParams(iss);

    //verify that they are equal
    ASSERT_EQ(params.level_expo, res_params->level_expo);
    ASSERT_EQ(params.plaintext_expo, res_params->plaintext_expo);
    ASSERT_EQ(params.torus_limbs, res_params->torus_limbs);
}

TEST(SERIALIZE_TEST, BigTorus) {
    BigTorusParams params(123, random(), random());
    BigTorus value(params);
    random(value);

    //serialize it into a string (instead of a file)
    ostringstream oss;
    serializeBigTorus(oss, value);
    string serial_result = oss.str();

    //deserialize
    istringstream iss(serial_result);
    shared_ptr<BigTorus> res = deserializeBigTorus(iss);

    //verify that they are equal
    ASSERT_LE(log2Diff(value, *res), -int64_t(value.params.torus_limbs * BITS_PER_LIMBS));
}

TEST(SERIALIZE_TEST, BigTorusVector) {
    BigTorusParams params(123, random(), random());
    int64_t length = 10;
    BigTorusVector value(length, params);

    for (int64_t i = 0; i < length; ++i) {
        random(value.getAT(i));
    }

    //serialize it into a string (instead of a file)
    ostringstream oss;
    serializeBigTorusVector(oss, value);
    string serial_result = oss.str();

    //deserialize
    istringstream iss(serial_result);
    shared_ptr<BigTorusVector> res = deserializeBigTorusVector(iss);

    //verify that they are equal
    for (int i = 0; i < length; ++i) {
        ASSERT_LE(log2Diff(value.getAT(i), res->getAT(i)), -int64_t(value.btp.torus_limbs * BITS_PER_LIMBS));
    }
}

TEST(SERIALIZE_TEST, TLweParams) {
    BigTorusParams params(123, random(), random());
    UINT64 N = 10;
    TLweParams tLweParams(N, params);


    //serialize it into a string (instead of a file)
    ostringstream oss;
    serializeTLweParams(oss, tLweParams);
    string serial_result = oss.str();

    //deserialize
    istringstream iss(serial_result);
    shared_ptr<TLweParams> res = deserializeTLweParams(iss);

    //verify that they are equal
    ASSERT_EQ(tLweParams.N, res->N);
    ASSERT_EQ(tLweParams.fixp_params.level_expo, res->fixp_params.level_expo);
    ASSERT_EQ(tLweParams.fixp_params.plaintext_expo, res->fixp_params.plaintext_expo);
    ASSERT_EQ(tLweParams.fixp_params.torus_limbs, res->fixp_params.torus_limbs);
}


TEST(SERIALIZE_TEST, TLweKey) {
    int64_t N = 10;
    TLweKey key(N);

    for (int64_t i = 0; i < N; i++) {
        int8_t coef = rand();
        key.key[i] = coef;
    }
    //serialize it into a string (instead of a file)
    ostringstream oss;
    serializeTLweKey(oss, key, N);
    string serial_result = oss.str();

    //deserialize
    istringstream iss(serial_result);
    shared_ptr<TLweKey> res = deserializeTLweKey(iss, N);

    //verify that they are equal
    for (int i = 0; i < N; ++i) {
        ASSERT_EQ(key.key[i], res->key[i]);
    }
}

TEST(SERIALIZE_TEST, TLwe) {
    BigTorusParams params(123, random(), random());
    int64_t N = 10;
    TLweParams tLweParams(N, params);
    TLwe value(tLweParams);

    //there are N+1 coeffs in a TLWE!
    for (int64_t i = 0; i < N + 1; ++i) {
        random(value.getAT(i));
    }

    //serialize it into a string (instead of a file)
    ostringstream oss;
    serializeTLwe(oss, value);
    string serial_result = oss.str();

    //deserialize
    istringstream iss(serial_result);
    shared_ptr<TLwe> res = deserializeTLwe(iss);

    //verify that they are equal
    for (int i = 0; i < N + 1; ++i) {
        ASSERT_LE(log2Diff(value.getAT(i), res->getAT(i)), -int64_t(value.btp.torus_limbs * BITS_PER_LIMBS));
    }
}

TEST(SERIALIZE_TEST, TRLweParams) {
    BigTorusParams params(123, random(), random());
    UINT64 N = 10;
    TRLweParams tRLweParams(N, params);


    //serialize it into a string (instead of a file)
    ostringstream oss;
    serializeTRLweParams(oss, tRLweParams);
    string serial_result = oss.str();

    //deserialize
    istringstream iss(serial_result);
    shared_ptr<TLweParams> res = deserializeTRLweParams(iss);

    //verify that they are equal
    ASSERT_EQ(tRLweParams.N, res->N);
    ASSERT_EQ(tRLweParams.fixp_params.level_expo, res->fixp_params.level_expo);
    ASSERT_EQ(tRLweParams.fixp_params.plaintext_expo, res->fixp_params.plaintext_expo);
    ASSERT_EQ(tRLweParams.fixp_params.torus_limbs, res->fixp_params.torus_limbs);
}

TEST(SERIALIZE_TEST, TRLwe) {
    BigTorusParams params(123, random(), random());
    int64_t N = 10;
    TRLweParams tLRweParams(N, params);
    TRLwe value(tLRweParams);

    //there are N+1 coeffs in a TLWE!
    for (int64_t i = 0; i < N; ++i) {
        random(value.a[0], 10);
        random(value.a[1], 10);
    }

    //serialize it into a string (instead of a file)
    ostringstream oss;
    serializeTRLwe(oss, value);
    string serial_result = oss.str();

    //deserialize
    istringstream iss(serial_result);
    shared_ptr<TRLwe> res = deserializeTRLwe(iss);

    //verify that they are equal
    for (int i = 0; i < N; ++i) {
        ASSERT_LE(log2Diff(value.a[0].getAT(i), res->a[0].getAT(i)),
                  -int64_t(value.params.fixp_params.torus_limbs * BITS_PER_LIMBS));
        ASSERT_LE(log2Diff(value.a[1].getAT(i), res->a[1].getAT(i)),
                  -int64_t(value.params.fixp_params.torus_limbs * BITS_PER_LIMBS));
    }
}

TEST(SERIALIZE_TEST, pubKsKey32) {

    BigTorusParams bt_params_in(123, random(), random());
    BigTorusParams bt_params_out(130, random(), random());


    int64_t N_in = 10;
    int64_t N_out = 8;
    int64_t alpha_bits = 80;

    TLweParams params_in(N_in, bt_params_in);
    TRLweParams params_out(N_out, bt_params_out);


    shared_ptr<TLweKey> key_in = tlwe_keygen(params_in);
    shared_ptr<TLweKey> key_out = tlwe_keygen(params_out);


    shared_ptr<pubKsKey32> ks_key = ks_keygen32(params_out, params_in, *key_in, *key_out, alpha_bits);


    //serialize it into a string (instead of a file)
    ostringstream oss;
    serializepubKsKey32(oss, *ks_key);
    string serial_result = oss.str();

    //deserialize
    istringstream iss(serial_result);
    shared_ptr<pubKsKey32> res = deserializepubKsKey32(iss);

    //verify that they are equal
    for (int64_t i = 0; i < int64_t(ks_key->in_params.N); i++) {
        for (int64_t j = 0; j < int64_t(ks_key->l_dec); j++) {
            for (int k = 0; k < ks_key->out_params.N; ++i) {
                ASSERT_LE(log2Diff(ks_key->kskey[i][j].a[0].getAT(k), res->kskey[i][j].a[0].getAT(k)),
                          -int64_t(ks_key->ks_params.fixp_params.torus_limbs * BITS_PER_LIMBS));

                ASSERT_LE(log2Diff(ks_key->kskey[i][j].a[1].getAT(k), res->kskey[i][j].a[1].getAT(k)),
                          -int64_t(ks_key->ks_params.fixp_params.torus_limbs * BITS_PER_LIMBS));

            }
        }
    }


}
