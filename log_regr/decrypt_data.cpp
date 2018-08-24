#include "io.h"
#include "common.h"
#include "lr_params.h"
#include "keyset.h"
#include "io_ctxt.h"

#include "tfhe_core.h"
#include "tfhe_gate.h"
#include "tfhe_functions.h"
#include "tfhe_alloc.h"


int main(int argc, char const *argv[]) {
    Data data;
    fill_data(data);

    LRParams lr_params;

    RandomGen::set_seed(42);

    const TfheParamSet *params = TfheParamSet::read(lr_params.params_filename);
    const TfheSecretKeySet *keyset = TfheSecretKeySet::read(lr_params.secret_keyset_filename, params);

    TGswSample<Torus>* X_cols_l1 = nullptr;
    TGswSample<Torus>* X_cols_l2 = nullptr;
    TLweSample<Torus>* y = nullptr;
    TLweSample<Torus>* sigmoid_xt_tps = nullptr;

    read_data(lr_params, sigmoid_xt_tps, y, X_cols_l1, X_cols_l2, params);

    {
        const TLweKey<Torus> *key = keyset->trlwe_key_l2;
        const TLweParams<Torus> *params = key->params;
        const int N = params->N;

        TorusPolynomial<Torus> *test_poly = new_obj<TorusPolynomial<Torus>>(N);

        for (int i = 0; i < 1; ++i) {
            TorusPolyFunctions<Torus>::Clear(test_poly);
            TLweFunctions<Torus>::Phase(test_poly, sigmoid_xt_tps+i, key);
            printf("sigmoid_xt_%d=\n[", i);
            for (int j = 0; j < N; ++j)
                printf("%lf, ", TorusUtils<Torus>::to_double(test_poly->coefsT[j]));
            printf("]\n");
        }
    }

    {
        const TLweKey<Torus> *key = keyset->trlwe_key_l2;
        const TLweParams<Torus> *params = key->params;
        const int N = params->N;

        TorusPolynomial<Torus> *msg = new_obj<TorusPolynomial<Torus>>(N);
        TorusPolyFunctions<Torus>::Clear(msg);

        TLweFunctions<Torus>::Phase(msg, y, key);

        printf("y:\n");
        for (int j = 0; j < data.n; ++j)
            printf("%lf ", TorusUtils<Torus>::to_double(msg->coefsT[j]));
        printf("\n");
    }

    {
        const TGswKey<Torus> *key = keyset->trgsw_key_l2;
        const TGswParams<Torus> *params = key->params;
        const int N = params->tlwe_params->N;

        IntPolynomial *msg = new_obj<IntPolynomial>(N);

        for (int j = 0; j < data.k; ++j) {
            TGswFunctions<Torus>::SymDecrypt(msg, X_cols_l1, key, 1<<lr_params.X_precbit);

            #ifndef NDEBUG
                printf("Xt[%d] (rescaled to int):\n", j);
                for (int i = 0; i < data.n; ++i)
                    printf("%d ", msg->coefs[i]);
                printf("\n");
            #endif
        }
    }

}
