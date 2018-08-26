#include "lr_params.h"
#include "tfhe_core.h"
#include "keyset.h"
#include "tlwe-functions-extra.h"
#include "io_ctxt.h"

using namespace std;

const int n_l0 = 666; // not on purpose, Nicola's formula was noise_alpha * 37 (here we have for 2^-18 noise std)
const int n_l1 = 2048;
const int n_l2 = 8192;
const int k = 1;

const double stddev_l0 = pow(2., -18);
const double stddev_l1 = pow(2., -48);
const double stddev_l2 = pow(2., -48);

const int l1_l = 8;
const int l1_Bgbit = 6;

const int l2_l = 8;
const int l2_Bgbit = 6;

const int ks_l1_l0_t = 3;
const int ks_l1_l0_basebit = 5;

const int ks_l2_l1_t = 5;
const int ks_l2_l1_basebit = 8;

TfheParamSet *new_parameters() {
    LweParams<Torus>* tlwe_params_l0 = new_obj<LweParams<Torus>>(n_l0, stddev_l0, -1);

    TLweParams<Torus>* trlwe_params_l1 = new_obj<TLweParams<Torus>>(n_l1, k, stddev_l1, -1);
    TGswParams<Torus>* trgsw_params_l1 = new_obj<TGswParams<Torus>>(l1_l, l1_Bgbit, trlwe_params_l1);

    TLweParams<Torus>* trlwe_params_l2 = new_obj<TLweParams<Torus>>(n_l2, k, stddev_l2, -1);
    TGswParams<Torus>* trgsw_params_l2 = new_obj<TGswParams<Torus>>(l2_l, l2_Bgbit, trlwe_params_l2);

    return new TfheParamSet(tlwe_params_l0, trgsw_params_l1, trgsw_params_l2);
}

TfheSecretKeySet* new_random_secret_keyset(const TfheParamSet *params) {
    LweKey<Torus>* tlwe_key_l0 = new_LweKey<Torus>(params->tlwe_params_l0);
    lweKeyGen(tlwe_key_l0);

    TGswKey<Torus>* trgsw_key_l1 = new_TGswKey<Torus>(params->trgsw_params_l1);
    tGswKeyGen(trgsw_key_l1);

    TGswKey<Torus>* trgsw_key_l2 = new_TGswKey<Torus>(params->trgsw_params_l2);
    tGswKeyGen(trgsw_key_l2);

    return new TfheSecretKeySet(params, tlwe_key_l0, trgsw_key_l1, trgsw_key_l2);
}

TGswSample<Torus>* create_bk(const LweKey<Torus>* tlwe_key, const TGswKey<Torus>* bk_key) {
    const int n = tlwe_key->params->n;
    const int* lwe_key_coefs = tlwe_key->key;
    const TGswParams<Torus>* bk_params = bk_key->params;

    double alpha = bk_params->tlwe_params->alpha_min;

    printf("Creating bk from l0 to l2:");
    printf("%8d input size", n);
    printf("%8d output size\n", bk_params->tlwe_params->N);

    TGswSample<Torus>* bk = new_obj_array<TGswSample<Torus>>(n, bk_params);

    #pragma omp parallel for
    for (int i=0; i<n; ++i) {
        TGswFunctions<Torus>::SymEncryptInt(bk+i, lwe_key_coefs[i], alpha, bk_key);
    }

    return bk;
}

LweKeySwitchKey<Torus>* create_ks_l1_l0(const LweKey<Torus>* out_key, const LweKey<Torus>* inp_key_extr) {
    const LweParams<Torus>* extract_params = inp_key_extr->params;
    const LweParams<Torus>* out_params = out_key->params;
    const int N = extract_params->n;

    printf("Creating ks from l1 to l0:");
    printf("%8d input size", N);
    printf("%8d output size\n", out_params->n);

    LweKeySwitchKey<Torus>* ks_l1_l0 = new_obj<LweKeySwitchKey<Torus>>(N, ks_l1_l0_t, ks_l1_l0_basebit, out_params);
    LweFunctions<Torus>::CreateKeySwitchKey(ks_l1_l0, inp_key_extr, out_key);

    return ks_l1_l0;
}

TLweKeySwitchKey<Torus>* create_ks_l2_l1(const TLweKey<Torus>* out_key, const LweKey<Torus>* inp_key_extr) {
    const LweParams<Torus>* extract_params = inp_key_extr->params;
    const TLweParams<Torus>* out_params = out_key->params;
    const int N = extract_params->n;

    printf("Creating ks from l2 to l1:");
    printf("%8d input size", N);
    printf("%8dx%8d output size\n", out_params->N, out_params->k);

    TLweKeySwitchKey<Torus>* ks_l2_l1 = new_obj<TLweKeySwitchKey<Torus>>(N, ks_l2_l1_t, ks_l2_l1_basebit, out_params);
    TLweFunctionsExtra<Torus>::CreateKeySwitchKey(ks_l2_l1, inp_key_extr, out_key);

    return ks_l2_l1;
}


TfheCloudKeySet* new_cloud_keyset(const TfheSecretKeySet* secret_keyset, const LRParams& lr_params) {
    TGswSample<Torus>* bk = create_bk(secret_keyset->tlwe_key_l0, secret_keyset->trgsw_key_l2);
    LweKeySwitchKey<Torus>* ks_l1_l0 = create_ks_l1_l0(secret_keyset->tlwe_key_l0, secret_keyset->tlwe_key_l1);
    TLweKeySwitchKey<Torus>* ks_l2_l1 = create_ks_l2_l1(secret_keyset->trlwe_key_l1, secret_keyset->tlwe_key_l2);;
    return new TfheCloudKeySet(secret_keyset->params, bk, ks_l1_l0, ks_l2_l1, false);
}


int main(int argc, char const *argv[]) {
    // uint32_t values[2];
    // values[0] = time(NULL)>>32;
    // values[1] = time(NULL);
    // RandomGen::set_seed(values, 2);
    RandomGen::set_seed(42);

    LRParams lr_params;

    // generate params
    TfheParamSet* params = new_parameters();
    TfheParamSet::write(lr_params.filename_params, params);

    // generate the secret keyset
    const TfheSecretKeySet *secret_keyset = new_random_secret_keyset(params);
    TfheSecretKeySet::write(lr_params.filename_secret_keyset, secret_keyset);

    // generate the cloud keyset
    const TfheCloudKeySet *cloud_keyset = new_cloud_keyset(secret_keyset, lr_params);
    TfheCloudKeySet::write(lr_params.filename_cloud_keyset, cloud_keyset);
}
