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

    const TfheParamSet *params = TfheParamSet::read(lr_params.params_filename);
    const TfheSecretKeySet *keyset = TfheSecretKeySet::read(lr_params.secret_keyset_filename, params);

    LweSample<Torus>* beta = read_tlwe_samples(argv[1], params->tlwe_params_l1, lr_params.k);
    for (int j = 0; j < lr_params.k; ++j) {
        Torus b = LweFunctions<Torus>::Phase(beta+j, keyset->tlwe_key_l1);
        printf("beta_%02d: %ld %lf\n", j, b, TorusUtils<Torus>::to_double(b));

        // static TORUS Phase(const LweSample<TORUS>* sample, const LweKey<TORUS>* key);
        // static TORUS SymDecrypt(const LweSample<TORUS>* sample, const LweKey<TORUS>* key, const int Msize);

    }
}
