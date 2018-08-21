#include "lr_params.h"
#include "tfhe_core.h"
#include "keyset.h"

using namespace std;

TfheParamSet *new_parameters() {
    const int n_l0 = 500;
    const int n_l1 = 2048;
    const int n_l2 = 8192;
    const int k = 1;

    const double stddev_l0 = pow(2., -15);
    const double stddev_l1 = pow(2., -48);
    const double stddev_l2 = pow(2., -48);

    const int l1_l = 8;
    const int l1_Bgbit = 6;

    const int l2_l = 8;
    const int l2_Bgbit = 6;


    LweParams<Torus>* tlwe_params_l0 = new_obj<LweParams<Torus>>(n_l0, stddev_l0, -1);

    TLweParams<Torus>* trlwe_params_l1 = new_obj<TLweParams<Torus>>(n_l1, k, stddev_l1, -1);
    TGswParams<Torus>* trgsw_params_l1 = new_obj<TGswParams<Torus>>(l1_l, l1_Bgbit, trlwe_params_l1);

    TLweParams<Torus>* trlwe_params_l2 = new_obj<TLweParams<Torus>>(n_l2, k, stddev_l2, -1);
    TGswParams<Torus>* trgsw_params_l2 = new_obj<TGswParams<Torus>>(l2_l, l2_Bgbit, trlwe_params_l2);

    return new TfheParamSet(tlwe_params_l0, trgsw_params_l1, trgsw_params_l2);
}

TfheSecretKeySet* new_random(const TfheParamSet *params) {
    LweKey<Torus>* tlwe_key_l0 = new_LweKey<Torus>(params->tlwe_params_l0);
    lweKeyGen(tlwe_key_l0);

    TGswKey<Torus>* trgsw_key_l1 = new_TGswKey<Torus>(params->trgsw_params_l1);
    tGswKeyGen(trgsw_key_l1);

    TGswKey<Torus>* trgsw_key_l2 = new_TGswKey<Torus>(params->trgsw_params_l2);
    tGswKeyGen(trgsw_key_l2);

    return new TfheSecretKeySet(params, tlwe_key_l0, trgsw_key_l1, trgsw_key_l2);
}


int main(int argc, char const *argv[]) {
    uint32_t values[2];
    values[0] = time(NULL)>>32;
    values[1] = time(NULL);
    RandomGen::set_seed(values, 2);

    LRParams lr_params;

    // generate params
    TfheParamSet* params = new_parameters();
    TfheParamSet::write(lr_params.params_filename, params);

    // generate the secret keyset
    const TfheSecretKeySet *secret_keyset = new_random(params);
    TfheSecretKeySet::write(lr_params.secret_keyset_filename, secret_keyset);

    // generate the cloud keyset
    const TfheCloudKeySet *cloud_keyset = new TfheCloudKeySet(secret_keyset);
    TfheCloudKeySet::write(lr_params.cloud_keyset_filename, cloud_keyset);
}
