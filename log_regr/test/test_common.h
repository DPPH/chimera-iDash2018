#include <include/gtest/gtest.h>

// #include "io.h"
#include "common.h"
#include "../keyset.h"
#include "../lr_params.h"

class TestCommon : public ::testing::Test {
public:
    void SetUp() {
        if (params != nullptr) return;

        params = TfheParamSet::read(lr_params.params_filename);
        secret_keyset = TfheSecretKeySet::read(lr_params.secret_keyset_filename, params);
        cloud_keyset = TfheCloudKeySet::read(lr_params.cloud_keyset_filename, params);

        read_data(lr_params, sigmoid_xt_tps, y, X_cols_l1, X_cols_l2, params);
    }

    void TearDown() {
    }

protected:
    // Data data;
    LRParams lr_params;
    const TfheParamSet* params = nullptr;
    const TfheSecretKeySet *secret_keyset = nullptr;
    const TfheCloudKeySet *cloud_keyset = nullptr;

    TGswSample<Torus>* X_cols_l1 = nullptr;
    TGswSample<Torus>* X_cols_l2 = nullptr;
    TLweSample<Torus>* y = nullptr;
    TLweSample<Torus>* sigmoid_xt_tps = nullptr;
};
