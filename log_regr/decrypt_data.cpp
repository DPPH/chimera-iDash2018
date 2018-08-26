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
    LRParams lr_params;

    const TfheParamSet *params = TfheParamSet::read(lr_params.filename_params);
    const TfheSecretKeySet *secret_keyset = TfheSecretKeySet::read(lr_params.filename_secret_keyset, params);

    TGswSample<Torus>* X_cols_l1 = nullptr;
    TGswSample<Torus>* X_cols_l2 = nullptr;
    TLweSample<Torus>* y = nullptr;
    TLweSample<Torus>* sigmoid_xt_tps = nullptr;

    read_data(lr_params, sigmoid_xt_tps, y, X_cols_l1, X_cols_l2, params);

    printf("Input data:\n");

    printf("X_cols_l1:\n");
    for (int j = 0; j < lr_params.k; ++j) {
        printf("col %1d:\n", j);
        print_X_col(X_cols_l1+j, secret_keyset->trgsw_key_l1, lr_params);
        // print_trgsw_sample(X_cols_l1+j, secret_keyset->trgsw_key_l1, lr_params.n, 1./lr_params.X_range);
        printf("\n");
    }

    // printf("X_cols_l2:\n");
    // for (int j = 0; j < lr_params.k; ++j) {
    //     printf("col %1d:\n", j);
        // print_X_col(X_cols_l2+j, secret_keyset->trgsw_key_l2, lr_params);
    //     printf("\n");
    // }

    printf("y[::-1]:\n");
    print_trlwe_sample(y, secret_keyset->trlwe_key_l2, lr_params.n, 1/lr_params.y_scale);
    printf("\n");

    printf("sigmoid_xt_tps:\n");
    for (int i = 0; i < 0; ++i) {
        printf("tp %3d:\n", i);
        print_trlwe_sample(sigmoid_xt_tps+i, secret_keyset->trlwe_key_l2, params->trlwe_params_l2->N, 1/lr_params.y_scale);
        printf("\n");
    }

}
