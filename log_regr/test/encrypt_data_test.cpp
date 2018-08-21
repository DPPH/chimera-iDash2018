#include <include/gtest/gtest.h>

#include "test_common.h"

class EncryptDataTest: public TestCommon {

};

TEST_F(EncryptDataTest, y_ctxt) {
    ifstream inp("y.ctxt", ifstream::binary);
    StdIstream inp_stream = to_Istream(inp);

    const TLweKey<Torus> *key = secret_keyset->trlwe_key_l2;
    const TLweParams<Torus> *params = key->params;
    const int N = params->N;

    TorusPolynomial<Torus> *msg = new_obj<TorusPolynomial<Torus>>(N);
    TorusPolyFunctions<Torus>::Clear(msg);
    TLweSample<Torus> *enc_y = new_obj<TLweSample<Torus>>(params);

    IOFunctions<Torus>::read_tLweSample(inp_stream, enc_y, params);
    TLweFunctions<Torus>::Phase(msg, enc_y, key);

    for (int j = 0; j < data.n; ++j) {
        double decr = TorusUtils<Torus>::to_double(msg->coefsT[j]);
        decr = decr > 0.5 ? decr-1 : decr;
        double real = data.y[data.n-j-1] * lr_params.scale;
        printf("%lf %lf\n", decr, real);
        ASSERT_NEAR(decr, real, 1e-6);
    }

    inp.close();
}

TEST_F(EncryptDataTest, X_cols) {

  // EXPECT_EQ(q0_.size(), 0);
}



