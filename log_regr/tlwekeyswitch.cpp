//
// Created by malika on 8/12/18.
//

#include "tlwekeyswitch.h"


template<typename TORUS>
TLweKeySwitchKey<TORUS>::TLweKeySwitchKey(
        const int n,
        const int t,
        const int basebit,
        const TLweParams<TORUS>* out_params,
        TLweSample<TORUS>* ks0_raw)
        :
        n(n),
        t(t),
        basebit(basebit),
        base(1<<basebit),
        out_params(out_params),
        ks0_raw(ks0_raw)
{
    ks = new TLweSample<TORUS>*[n];
    for (int p = 0; p < n; ++p)
        ks[p] = ks0_raw + t*p;
}

template<typename TORUS>
TLweKeySwitchKey<TORUS>::~TLweKeySwitchKey() {
    delete[] ks;
}

/**
 * TLweKeySwitchKey constructor function
 */
template<typename TORUS>
void TLweKeySwitchKey<TORUS>::init_obj(TLweKeySwitchKey<TORUS>* obj, int n, int t, int basebit, const TLweParams<TORUS>* out_params) {
    TLweSample<TORUS>* ks0_raw = new_obj_array<TLweSample<TORUS>>(n*t, out_params);
    new(obj) TLweKeySwitchKey<TORUS>(n,t,basebit,out_params, ks0_raw);
}

/**
 * TLweKeySwitchKey destructor
 */
template<typename TORUS>
void TLweKeySwitchKey<TORUS>::destroy_obj(TLweKeySwitchKey<TORUS>* obj) {
    const int n = obj->n;
    const int t = obj->t;
    del_obj_array(n*t, obj->ks0_raw);
    obj->~TLweKeySwitchKey();
}


template<typename TORUS>
void TLweKeySwitchKey<TORUS>::write(Ostream& out_stream) const {
    out_stream.fwrite(&n, sizeof(int32_t));
    out_stream.fwrite(&t, sizeof(int32_t));
    out_stream.fwrite(&basebit, sizeof(int32_t));
    for (int i = 0; i < n*t; ++i) {
        IOFunctions<TORUS>::write_tLweSample(out_stream, ks0_raw+i, out_params);
    }
}

template<typename TORUS>
TLweKeySwitchKey<TORUS>* TLweKeySwitchKey<TORUS>::read_new(Istream& inp_stream, const TLweParams<TORUS>* params) {
    int32_t n, t, basebit;

    inp_stream.fread(&n, sizeof(int32_t));
    inp_stream.fread(&t, sizeof(int32_t));
    inp_stream.fread(&basebit, sizeof(int32_t));

    TLweKeySwitchKey<TORUS>* ks = new_obj<TLweKeySwitchKey<TORUS>>(n, t, basebit, params);
    for (int i = 0; i < n*t; ++i) {
        IOFunctions<TORUS>::read_tLweSample(inp_stream, ks->ks0_raw+i, params);
    }

    return ks;
}
