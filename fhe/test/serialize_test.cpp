#include <include/gtest/gtest.h>
#include "../BigTorus.h"
#include "../BigTorusVector.h"
#include "../TLwe.h"
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
