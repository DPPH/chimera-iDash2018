#ifndef KEYSET_H
#define KEYSET_H

#include "tfhe_core.h"
#include "tfhe_io.h"
#include "tgsw.h"

#define Torus Torus64

struct TfheParamSet {
    TfheParamSet(
        LweParams<Torus>* tlwe_params_l0,
        TGswParams<Torus>* trgsw_params_l1,
        TGswParams<Torus>* trgsw_params_l2
        )
    :
        tlwe_params_l0(tlwe_params_l0),
        trgsw_params_l1(trgsw_params_l1),
        trlwe_params_l1(trgsw_params_l1->tlwe_params),
        trgsw_params_l2(trgsw_params_l2),
        trlwe_params_l2(trgsw_params_l2->tlwe_params)
    { }

    static void write(Ostream& out, const TfheParamSet*);
    static void write(const char* filename, const TfheParamSet*);
    static const TfheParamSet* read(Istream& inp);
    static const TfheParamSet* read(const char* filename);

    const LweParams<Torus>* const tlwe_params_l0 = nullptr;

    const TGswParams<Torus>* const trgsw_params_l1 = nullptr;
    const TLweParams<Torus>* const trlwe_params_l1 = nullptr;   // reference to trgsw_params_l2->tlwe_params

    const TGswParams<Torus>* const trgsw_params_l2 = nullptr;
    const TLweParams<Torus>* const trlwe_params_l2 = nullptr;   // reference to trgsw_params_l2->tlwe_params

    double alpha_l0 = pow(2., -12);
    double alpha_l1 = pow(2., -48);
    double alpha_l2 = pow(2., -48);

    double alpha_test_poly;
    double alpha_y;
    double alpha_X_cols;

    double alpha_bk;

};


class TfheSecretKeySet {
public:
    TfheSecretKeySet(
        const TfheParamSet* params,
        const LweKey<Torus>* tlwe_key_l0,
        const TGswKey<Torus>* trgsw_key_l1,
        const TGswKey<Torus>* trgsw_key_l2
        )
    :
        params(params),
        tlwe_key_l0(tlwe_key_l0),
        trgsw_key_l1(trgsw_key_l1),
        trlwe_key_l1(&(trgsw_key_l1->tlwe_key)),
        trgsw_key_l2(trgsw_key_l2),
        trlwe_key_l2(&(trgsw_key_l2->tlwe_key))
    { }

    static void del(const TfheSecretKeySet*);

    static void write(Ostream& out, const TfheSecretKeySet* const);
    static void write(const char* filename, const TfheSecretKeySet* const);
    static const TfheSecretKeySet* read(Istream& inp, const TfheParamSet* params);
    static const TfheSecretKeySet* read(const char* filename, const TfheParamSet* params);


    const TfheParamSet* params;

    // level 0
    const LweKey<Torus>* tlwe_key_l0;

    // level 1
    const TGswKey<Torus>* trgsw_key_l1;
    const TLweKey<Torus>* trlwe_key_l1;  // reference to trgsw_key_l2->tlwe_key

    // level 2
    const TGswKey<Torus>* trgsw_key_l2;
    const TLweKey<Torus>* trlwe_key_l2; // reference to trgsw_key_l2->tlwe_key

};

class TfheCloudKeySet {
    void create_bk(const TfheSecretKeySet* secret_keyset);

    void create_ks_l1_l0(const TfheSecretKeySet* secret_keyset);

    void create_bk_fft();


public:
    TfheCloudKeySet(const TfheSecretKeySet* secret_keyset)
    :
        params(secret_keyset->params)
    {
        create_bk(secret_keyset);
        create_bk_fft();
        create_ks_l1_l0(secret_keyset);
    }

    TfheCloudKeySet(
        const TfheParamSet* params,
        TGswSample<Torus>* bk, //TRGSW encryption of sk_l0
        LweKeySwitchKey<Torus>* ks_l1_l0
        // const TLweKeySwitchKey<Torus>* ks_l2_l1,
        )
    :
        params(params),
        bk(bk),
        ks_l1_l0(ks_l1_l0)
        // ks_l2_l1(ks_l2_l1)
    {
        create_bk_fft();
    }

    static void del(const TfheCloudKeySet*);

    static void write(Ostream& out, const TfheCloudKeySet*);
    static void write(const char* filename, const TfheCloudKeySet*);
    static const TfheCloudKeySet* read(Istream& inp, const TfheParamSet* params);
    // static const TfheCloudKeySet* read(const char* filename, const TfheParamSet* params);


    const TfheParamSet* params;
    TGswSample<Torus>* bk = nullptr;
    TGswSampleFFT<Torus>* bk_fft = nullptr;
    LweKeySwitchKey<Torus>* ks_l1_l0 = nullptr;
    // const TLweKeySwitchKey<Torus>* ks_l2_l1 = nullptr;
};

#endif

