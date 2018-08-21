#include <include/gtest/gtest.h>

// #include "io.h"
#include "common.h"
#include "../keyset.h"
#include "../lr_params.h"

class TestCommon : public ::testing::Test {
public:
    void SetUp() {
        if (params != nullptr) return;

        // fill_data(data);

        // lr_params.n = data.n;
        // lr_params.k = data.k;
        // lr_params.m = data.m;

        params = TfheParamSet::read(lr_params.params_filename);
        secret_keyset = TfheSecretKeySet::read(lr_params.secret_keyset_filename, params);
        cloud_keyset = TfheCloudKeySet::read(lr_params.cloud_keyset_filename, params);
    }

    void TearDown() {
    }

protected:
    // Data data;
    LRParams lr_params;
    const TfheParamSet* params = nullptr;
    const TfheSecretKeySet *secret_keyset = nullptr;
    const TfheCloudKeySet *cloud_keyset = nullptr;
};
