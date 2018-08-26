#include "keyset.h"

#include "tfhe_core.h"
#include "tfhe_io.h"

#include <cassert>
#include <iostream>
#include <fstream>

using namespace std;

void TfheParamSet::write(Ostream& out, const TfheParamSet* params) {
    IOFunctions<Torus>::write_lweParams(out, params->tlwe_params_l0);
    IOFunctions<Torus>::write_tGswParams(out, params->trgsw_params_l1);
    IOFunctions<Torus>::write_tGswParams(out, params->trgsw_params_l2);
}

void TfheParamSet::write(const char* filename, const TfheParamSet* params) {
    ofstream out(filename, ofstream::binary);
    if (out.is_open()) {
        StdOstream out_stream = to_Ostream(out);
        TfheParamSet::write(out_stream, params);
        out.close();
    } else {
        fprintf(stderr, "Cannot open file '%s' for writing tfhe parameters\n", filename);
    }
}

const TfheParamSet* TfheParamSet::read(Istream& inp) {
    LweParams<Torus>* tlwe_params_l0 = IOFunctions<Torus>::read_new_lweParams(inp);
    TGswParams<Torus>* trgsw_params_l1 = IOFunctions<Torus>::read_new_tGswParams(inp);
    TGswParams<Torus>* trgsw_params_l2 = IOFunctions<Torus>::read_new_tGswParams(inp);

    return new TfheParamSet(tlwe_params_l0, trgsw_params_l1, trgsw_params_l2);
}

const TfheParamSet* TfheParamSet::read(const char* filename) {
    ifstream inp(filename, ifstream::binary);
    const TfheParamSet* params = nullptr;
    if (inp.is_open()) {
        StdIstream inp_stream = to_Istream(inp);
        params = TfheParamSet::read(inp_stream);
        inp.close();
    } else {
        fprintf(stderr, "Cannot open parameters file '%s'\n", filename);
    }
    return params;
}

void TfheCloudKeySet::create_bk_fft() {
    assert(params != nullptr);
    assert(bk != nullptr);

    const int n = params->tlwe_params_l0->n;
    const TGswParams<Torus>* bk_params = params->trgsw_params_l2;

    // Bootstrapping Key FFT
    bk_fft = new_obj_array<TGswSampleFFT<Torus>>(n, bk_params);

    #ifndef NDEBUG
        printf("Creating bk_fft from l0 to l2\n");
    #endif

    // #pragma omp parallel for // tGswToFFTConvert is not thread safe?!?
    for (int i=0; i<n; ++i) {
        TGswFunctions<Torus>::ToFFTConvert(bk_fft+i, bk+i, bk_params);
    }
}


void TfheSecretKeySet::write(Ostream& out, const TfheSecretKeySet* secret_keyset) {
    IOFunctions<Torus>::write_lweKey(out, secret_keyset->tlwe_key_l0, false);
    IOFunctions<Torus>::write_tGswKey(out, secret_keyset->trgsw_key_l1, false);
    IOFunctions<Torus>::write_tGswKey(out, secret_keyset->trgsw_key_l2, false);
}

void TfheSecretKeySet::write(const char* filename, const TfheSecretKeySet* secret_keyset) {
    ofstream out(filename, ofstream::binary);
    if (out.is_open()) {
        StdOstream out_stream = to_Ostream(out);
        TfheSecretKeySet::write(out_stream, secret_keyset);
        out.close();
    } else {
        fprintf(stderr, "Cannot open file '%s' for writing secret keyset\n", filename);
    }
}


const TfheSecretKeySet* TfheSecretKeySet::read(Istream& inp, const TfheParamSet* params) {
    LweKey<Torus>* tlwe_key_l0 = IOFunctions<Torus>::read_new_lweKey(inp, params->tlwe_params_l0);
    TGswKey<Torus>* trgsw_key_l1 = IOFunctions<Torus>::read_new_tGswKey(inp, params->trgsw_params_l1);
    TGswKey<Torus>* trgsw_key_l2 = IOFunctions<Torus>::read_new_tGswKey(inp, params->trgsw_params_l2);

    return new TfheSecretKeySet(params, tlwe_key_l0, trgsw_key_l1, trgsw_key_l2);
}

const TfheSecretKeySet* TfheSecretKeySet::read(const char* filename, const TfheParamSet* params) {
    ifstream inp(filename, ifstream::binary);
    const TfheSecretKeySet* keyset = nullptr;
    if (inp.is_open()) {
        StdIstream inp_stream = to_Istream(inp);
        keyset = TfheSecretKeySet::read(inp_stream, params);
        inp.close();
    } else {
        fprintf(stderr, "Cannot open secret keyset file '%s'\n", filename);
    }
    return keyset;
}

void TfheCloudKeySet::write(Ostream& out, const TfheCloudKeySet* keyset) {
    const TfheParamSet* params = keyset->params;

    const int n = params->tlwe_params_l0->n;
    const TGswParams<Torus>* bk_params = params->trgsw_params_l2;
    const TGswSample<Torus>* bk = keyset->bk;
    for (int i = 0; i < n; ++i)
        IOFunctions<Torus>::write_tGswSample(out, bk+i, bk_params);

    IOFunctions<Torus>::write_lweKeySwitchKey(out, keyset->ks_l1_l0, false);
    keyset->ks_l2_l1->write(out);
}

void TfheCloudKeySet::write(const char* filename, const TfheCloudKeySet* secret_keyset) {
    ofstream out(filename, ofstream::binary);
    if (out.is_open()) {
        StdOstream out_stream = to_Ostream(out);
        TfheCloudKeySet::write(out_stream, secret_keyset);
        out.close();
    } else {
        fprintf(stderr, "Cannot open file '%s' for writing cloud keyset\n", filename);
    }
}

const TfheCloudKeySet* TfheCloudKeySet::read(Istream& inp, const TfheParamSet* params) {
    const int n = params->tlwe_params_l0->n;
    const TGswParams<Torus>* bk_params = params->trgsw_params_l2;
    TGswSample<Torus>* bk = new_obj_array<TGswSample<Torus>>(n, bk_params);
    for (int i = 0; i < n; ++i)
        IOFunctions<Torus>::read_tGswSample(inp, bk+i, bk_params);

    LweKeySwitchKey<Torus>* ks_l1_l0 = IOFunctions<Torus>::read_new_lweKeySwitchKey(inp, params->tlwe_params_l0);
    TLweKeySwitchKey<Torus>* ks_l2_l1 = TLweKeySwitchKey<Torus>::read_new(inp, params->trlwe_params_l1);

    return new TfheCloudKeySet(params, bk, ks_l1_l0, ks_l2_l1);
}

const TfheCloudKeySet* TfheCloudKeySet::read(const char* filename, const TfheParamSet* params) {
    ifstream inp(filename, ifstream::binary);
    const TfheCloudKeySet* keyset = nullptr;
    if (inp.is_open()) {
        StdIstream inp_stream = to_Istream(inp);
        keyset = TfheCloudKeySet::read(inp_stream, params);
        inp.close();
    } else {
        fprintf(stderr, "Cannot open cloud keyset file '%s'\n", filename);
    }
    return keyset;
}
