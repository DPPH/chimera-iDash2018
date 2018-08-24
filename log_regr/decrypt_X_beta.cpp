#include "lr_params.h"
#include "keyset.h"
#include "io_ctxt.h"
#include "io.h"

#include "tfhe_core.h"
#include "tfhe_gate.h"
#include "tfhe_functions.h"
#include "tfhe_alloc.h"

int main(int argc, char const *argv[]) {
    LRParams lr_params;
    read_data_header(lr_params);
    lr_params.update();

    const TfheParamSet *params = TfheParamSet::read(lr_params.params_filename);
    const TfheSecretKeySet *secret_keyset = TfheSecretKeySet::read(lr_params.secret_keyset_filename, params);

    LweSample<Torus>* X_beta = read_tlwe_samples(argv[1], params->tlwe_params_l0, lr_params.n);

    for (int i = 0; i < lr_params.n; ++i) {
        print_tlwe_sample(X_beta+i, secret_keyset->tlwe_key_l0, 1./lr_params.X_beta_scale);
    }
    printf("\n");
}
