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

void encrypt_lwe(const char *filename, const LweKey<Torus> *key) {
    const LweParams<Torus> *params = key->params;
    // const double alpha = params->alpha_min;
    const double alpha = 0.0;

    LweSample<Torus> *msg_enc = new_obj<LweSample<Torus>>(params);
    printf("enc %s: %ld\n", filename, 171717171L);

    LweFunctions<Torus>::SymEncrypt(msg_enc, 171717171L, alpha, key);

    ofstream out(filename, ofstream::binary);
    StdOstream out_stream = to_Ostream(out);
    IOFunctions<Torus>::write_lweSample(out_stream, msg_enc, params);
    out.close();
}

void decrypt_lwe(const char *filename, const LweKey<Torus> *key) {
    const LweParams<Torus> *params = key->params;

    LweSample<Torus> *msg_enc = new_obj<LweSample<Torus>>(params);

    ifstream inp(filename, ifstream::binary);
    StdIstream inp_stream = to_Istream(inp);
    IOFunctions<Torus>::read_lweSample(inp_stream, msg_enc, params);
    inp.close();

    printf("dec %s: %ld\n", filename, LweFunctions<Torus>::Phase(msg_enc, key));
}

void encrypt_tlwe(const char *filename, const TLweKey<Torus> *key) {
    const TLweParams<Torus> *params = key->params;
    const int N = params->N;
    // const double alpha = params->alpha_min;
    const double alpha = 0.0;

    TorusPolynomial<Torus> *msg = new_obj<TorusPolynomial<Torus>>(N);
    TLweSample<Torus> *msg_enc = new_obj<TLweSample<Torus>>(params);

    TorusPolyFunctions<Torus>::Clear(msg);
    msg->coefsT[0] = 1300001;
    msg->coefsT[1] = 171717171;

    printf("enc %s: %ld %ld\n", filename, msg->coefsT[0], msg->coefsT[1]);
    TLweFunctions<Torus>::SymEncrypt(msg_enc, msg, alpha, key);

    ofstream out(filename, ofstream::binary);
    StdOstream out_stream = to_Ostream(out);
    IOFunctions<Torus>::write_tLweSample(out_stream, msg_enc, params);
    out.close();
}

void decrypt_tlwe(const char *filename, const TLweKey<Torus> *key) {
    const TLweParams<Torus> *params = key->params;
    const int N = params->N;

    TorusPolynomial<Torus> *msg = new_obj<TorusPolynomial<Torus>>(N);
    TLweSample<Torus> *msg_enc = new_obj<TLweSample<Torus>>(params);
    TorusPolyFunctions<Torus>::Clear(msg);

    ifstream inp(filename, ifstream::binary);
    StdIstream inp_stream = to_Istream(inp);
    IOFunctions<Torus>::read_tLweSample(inp_stream, msg_enc, params);
    inp.close();

    TLweFunctions<Torus>::Phase(msg, msg_enc, key);

    printf("dec %s: %ld %ld\n", filename, msg->coefsT[0], msg->coefsT[1]);
}

void encrypt_tgsw(const char *filename, const TGswKey<Torus> *key) {
    const TGswParams<Torus> *params = key->params;
    const int N = params->tlwe_params->N;
    // const double alpha = params->alpha_min;
    const double alpha = 0.0;

    IntPolynomial *msg = new_obj<IntPolynomial>(N);
    TGswSample<Torus>* msg_enc = new_obj<TGswSample<Torus>>(params);

    IntPolyFunctions::Clear(msg);
    msg->coefs[0] = 1300001;
    msg->coefs[1] = 171717171;

    printf("enc %s: %d %d\n", filename, msg->coefs[0], msg->coefs[1]);
    TGswFunctions<Torus>::SymEncrypt(msg_enc, msg, alpha, key);

    ofstream out(filename, ofstream::binary);
    StdOstream out_stream = to_Ostream(out);
    IOFunctions<Torus>::write_tGswSample(out_stream, msg_enc, params);
    out.close();
}

void decrypt_tgsw(const char *filename, const TGswKey<Torus> *key) {
    const TGswParams<Torus> *params = key->params;
    const int N = params->tlwe_params->N;

    IntPolynomial *msg = new_obj<IntPolynomial>(N);
    TGswSample<Torus>* msg_enc = new_obj<TGswSample<Torus>>(params);

    IntPolyFunctions::Clear(msg);

    ifstream inp(filename, ifstream::binary);
    StdIstream inp_stream = to_Istream(inp);
    IOFunctions<Torus>::read_tGswSample(inp_stream, msg_enc, params);
    inp.close();

    TGswFunctions<Torus>::SymDecrypt(msg, msg_enc, key, 1L<<30);

    printf("enc %s: %d %d\n", filename, msg->coefs[0], msg->coefs[1]);
}

int main() {
    Data data;
    fill_data(data);

    LRParams lr_params;
    lr_params.n = data.n;
    lr_params.k = data.k;
    lr_params.m = data.m;

    RandomGen::set_seed(42);

    const TfheParamSet *params = TfheParamSet::read(lr_params.params_filename);
    const TfheSecretKeySet *keyset = TfheSecretKeySet::read(lr_params.secret_keyset_filename, params);


    encrypt_lwe("tlwe.l0", keyset->tlwe_key_l0);
    decrypt_lwe("tlwe.l0", keyset->tlwe_key_l0);

    encrypt_tgsw("trgsw.l2", keyset->trgsw_key_l2);
    decrypt_tgsw("trgsw.l2", keyset->trgsw_key_l2);

    encrypt_tgsw("trgsw.l1", keyset->trgsw_key_l1);
    decrypt_tgsw("trgsw.l1", keyset->trgsw_key_l1);

    encrypt_tlwe("trlwe.l1", keyset->trlwe_key_l1);
    decrypt_tlwe("trlwe.l1", keyset->trlwe_key_l1);

    encrypt_tlwe("trlwe.l2", keyset->trlwe_key_l2);
    decrypt_tlwe("trlwe.l2", keyset->trlwe_key_l2);


    return 0;
}
