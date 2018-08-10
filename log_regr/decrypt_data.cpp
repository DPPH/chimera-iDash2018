#include "io.h"
#include "common.h"
#include "params.h"
#include "keyset.h"

#include "tfhe_core.h"
#include "tfhe_gate.h"
#include "tfhe_functions.h"
#include "tfhe_alloc.h"
#include "tfhe_io.h"

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char const *argv[]) {
    Data data;
    fill_data(data);

    LRParams lr_params;
    lr_params.n = data.n;
    lr_params.k = data.k;
    lr_params.m = data.m;

    RandomGen::set_seed(42);

    const TfheParamSet *params = TfheParamSet::read(lr_params.params_filename);
    const TfheSecretKeySet *keyset = TfheSecretKeySet::read(lr_params.secret_keyset_filename, params);

    {
        ifstream inp("sigmoid_xt.test_poly", ifstream::binary);
        StdIstream inp_stream = to_Istream(inp);

        const TLweKey<Torus> *key = keyset->trlwe_key_l2;
        const TLweParams<Torus> *params = key->params;
        const int N = params->N;

        TorusPolynomial<Torus> *test_poly = new_obj<TorusPolynomial<Torus>>(N);
        TLweSample<Torus> *test_poly_enc = new_obj<TLweSample<Torus>>(params);

        for (int i = 0; i < 1; ++i) {
            IOFunctions<Torus>::read_tLweSample(inp_stream, test_poly_enc, params);
            TorusPolyFunctions<Torus>::Clear(test_poly);
            TLweFunctions<Torus>::Phase(test_poly, test_poly_enc, key);
            printf("sigmoid_xt[%d]:\n", i);
            for (int j = 0; j < N; ++j)
                printf("%lf ", TorusUtils<Torus>::to_double(test_poly->coefsT[j]));
            printf("\n");
        }

        inp.close();
    }

    {
        ifstream inp("y.ctxt", ifstream::binary);
        StdIstream inp_stream = to_Istream(inp);

        const TLweKey<Torus> *key = keyset->trlwe_key_l1;
        const TLweParams<Torus> *params = key->params;
        const int N = params->N;

        TorusPolynomial<Torus> *msg = new_obj<TorusPolynomial<Torus>>(N);
        TorusPolyFunctions<Torus>::Clear(msg);
        TLweSample<Torus> *enc_y = new_obj<TLweSample<Torus>>(params);

        IOFunctions<Torus>::read_tLweSample(inp_stream, enc_y, params);
        TLweFunctions<Torus>::Phase(msg, enc_y, key);

        printf("y:\n");
        for (int j = 0; j < data.n; ++j)
            printf("%lf ", TorusUtils<Torus>::to_double(msg->coefsT[j]));
        printf("\n");

        inp.close();
    }

    {
        ifstream inp("X_cols.ctxt", ifstream::binary);
        StdIstream inp_stream = to_Istream(inp);

        const TGswKey<Torus> *key = keyset->trgsw_key_l1;
        const TGswParams<Torus> *params = key->params;
        const int N = params->tlwe_params->N;

        IntPolynomial *msg = new_obj<IntPolynomial>(N);
        TGswSample<Torus>* enc_X_col = new_obj<TGswSample<Torus>>(params);

        for (int j = 0; j < data.k; ++j) {
            IOFunctions<Torus>::read_tGswSample(inp_stream, enc_X_col, params);
            TGswFunctions<Torus>::SymDecrypt(msg, enc_X_col, key, 1<<lr_params.X_precbit);

            #ifndef NDEBUG
                printf("Xt[%d] (rescaled to int):\n", j);
                for (int i = 0; i < data.n; ++i)
                    printf("%d ", msg->coefs[i]);
                printf("\n");
            #endif
        }

        inp.close();
    }

}
