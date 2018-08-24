#include "lr_params.h"
#include "keyset.h"
#include "io_ctxt.h"

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

    LweSample<Torus>* beta = read_tlwe_samples(argv[1], params->tlwe_params_l2, lr_params.k);
    double scale = 1/lr_params.beta_scale;
    // printf("beta: ");
    for (int j = 0; j < lr_params.k; ++j) {
        print_tlwe_sample(beta+j, secret_keyset->tlwe_key_l2, scale);
    }
    printf("\n");
}
