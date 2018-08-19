#include <include/gtest/gtest.h>
#include "../BigTorus.h"
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
