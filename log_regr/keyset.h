#ifndef KEYSET_H
#define KEYSET_H

#include "tfhe_core.h"
#include "tfhe_io.h"
// #include "tgsw.h"
#include "tlwekeyswitch.h"

#include <cassert>

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
        tlwe_params_l1(&(trlwe_params_l1->extracted_lweparams)),
        trgsw_params_l2(trgsw_params_l2),
        trlwe_params_l2(trgsw_params_l2->tlwe_params),
        tlwe_params_l2(&(trlwe_params_l2->extracted_lweparams))
    { }

    ~TfheParamSet() {
        del_obj(tlwe_params_l0);
        del_obj(trgsw_params_l1);
        del_obj(trgsw_params_l2);
    }

    static void write(Ostream& out, const TfheParamSet*);
    static void write(const char* filename, const TfheParamSet*);
    static const TfheParamSet* read(Istream& inp);
    static const TfheParamSet* read(const char* filename);

    LweParams<Torus>* const tlwe_params_l0 = nullptr;

    TGswParams<Torus>* const trgsw_params_l1 = nullptr;
    const TLweParams<Torus>* const trlwe_params_l1 = nullptr;    // reference to trgsw_params_l1->tlwe_params
    const LweParams<Torus>* const tlwe_params_l1 = nullptr;      // extracted LWE params from TLWE L1

    TGswParams<Torus>* const trgsw_params_l2 = nullptr;
    const TLweParams<Torus>* const trlwe_params_l2 = nullptr;    // reference to trgsw_params_l2->tlwe_params
    const LweParams<Torus>* const tlwe_params_l2 = nullptr;      // extracted LWE params from TLWE L2
};


class TfheSecretKeySet {
public:
    TfheSecretKeySet(
        const TfheParamSet* params,
        LweKey<Torus>* tlwe_key_l0,
        TGswKey<Torus>* trgsw_key_l1,
        TGswKey<Torus>* trgsw_key_l2
        )
    :
        params(params),
        tlwe_key_l0(tlwe_key_l0),
        trgsw_key_l1(trgsw_key_l1),
        trlwe_key_l1(&(trgsw_key_l1->tlwe_key)),
        trgsw_key_l2(trgsw_key_l2),
        trlwe_key_l2(&(trgsw_key_l2->tlwe_key))
    {
        assert(tlwe_key_l0->params == params->tlwe_params_l0);
        assert(trgsw_key_l1->params == params->trgsw_params_l1);
        assert(trgsw_key_l2->params == params->trgsw_params_l2);

        tlwe_key_l1 = new_obj<LweKey<Torus>>(params->tlwe_params_l1);
        TLweFunctions<Torus>::ExtractKey(tlwe_key_l1, trlwe_key_l1);

        tlwe_key_l2 = new_obj<LweKey<Torus>>(params->tlwe_params_l2);
        TLweFunctions<Torus>::ExtractKey(tlwe_key_l2, trlwe_key_l2);
    }

    ~TfheSecretKeySet() {
        del_obj(tlwe_key_l2);
        del_obj(tlwe_key_l1);
        del_obj(trgsw_key_l2);
        del_obj(trgsw_key_l1);
        del_obj(tlwe_key_l0);
    }

    static void del(const TfheSecretKeySet*);

    static void write(Ostream& out, const TfheSecretKeySet* const);
    static void write(const char* filename, const TfheSecretKeySet* const);
    static const TfheSecretKeySet* read(Istream& inp, const TfheParamSet* params);
    static const TfheSecretKeySet* read(const char* filename, const TfheParamSet* params);


    const TfheParamSet* params = nullptr;

    // level 0
    LweKey<Torus>* tlwe_key_l0 = nullptr;

    // level 1
    TGswKey<Torus>* trgsw_key_l1 = nullptr;
    TLweKey<Torus>* trlwe_key_l1 = nullptr;  // reference to trgsw_key_l2->tlwe_key
    LweKey<Torus>* tlwe_key_l1 = nullptr;  // extracted TLWE key from trlwe_key_l1

    // level 2
    TGswKey<Torus>* trgsw_key_l2 = nullptr;
    TLweKey<Torus>* trlwe_key_l2 = nullptr; // reference to trgsw_key_l2->tlwe_key
    LweKey<Torus>* tlwe_key_l2 = nullptr;  // extracted TLWE key from trlwe_key_l2
};

class TfheCloudKeySet {
private:
    void create_bk_fft();

public:
    TfheCloudKeySet(
        const TfheParamSet* params,
        TGswSample<Torus>* bk, //TRGSW encryption of TLWE L0 SK
        LweKeySwitchKey<Torus>* ks_l1_l0, //KS of TLWE L1 SK to TLWE L0 SK
        TLweKeySwitchKey<Torus>* ks_l2_l1, //KS of TLWE L2 SK to TRLWE L1 SK
        bool compute_bk_fft = true
        )
    :
        params(params),
        bk(bk),
        ks_l1_l0(ks_l1_l0),
        ks_l2_l1(ks_l2_l1)
    {
        if (compute_bk_fft)
            create_bk_fft();
    }

    ~TfheCloudKeySet() {
        const int n = params->tlwe_params_l0->n;
        if (bk_fft != nullptr) del_obj_array(n, bk_fft);
        del_obj_array(n, bk);
        del_obj(ks_l2_l1);
        del_obj(ks_l1_l0);
    }


    static void del(const TfheCloudKeySet*);

    static void write(Ostream& out, const TfheCloudKeySet*);
    static void write(const char* filename, const TfheCloudKeySet*);
    static const TfheCloudKeySet* read(Istream& inp, const TfheParamSet* params);
    static const TfheCloudKeySet* read(const char* filename, const TfheParamSet* params);


    const TfheParamSet* const params = nullptr;
    TGswSample<Torus>* const bk = nullptr;
    LweKeySwitchKey<Torus>* const ks_l1_l0 = nullptr;
    TLweKeySwitchKey<Torus>* const ks_l2_l1 = nullptr;

    TGswSampleFFT<Torus>* bk_fft = nullptr;
};

#endif

