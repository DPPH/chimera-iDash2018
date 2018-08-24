/
// Created by malika on 8/12/18.
//

#include "tlwekeyswitch.h"
#include "tfhe_core.h"
#include "tfhe_alloc.h"
#include "tlwe.h"

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
    ks1_raw = new TLweSample<TORUS>*[n*t];
    ks = new TLweSample<TORUS>**[n];

    for (int p = 0; p < n*t; ++p)
        ks1_raw[p] = ks0_raw + base*p;
    for (int p = 0; p < n; ++p)
        ks[p] = ks1_raw + t*p;
}

template<typename TORUS>
TLweKeySwitchKey<TORUS>::~TLweKeySwitchKey() {
    delete[] ks;
    delete[] ks1_raw;
}

/**
 * TLweKeySwitchKey constructor function
 */
template<typename TORUS>
void TLweKeySwitchKey<TORUS>::init_obj(TLweKeySwitchKey<TORUS>* obj, int n, int t, int basebit, const TLweParams<TORUS>* out_params) {
    const int base=1<<basebit;
    TLweSample<TORUS>* ks0_raw = new_obj_array<TLweSample<TORUS>>(n*t*base, out_params);
    new(obj) TLweKeySwitchKey<TORUS>(n,t,basebit,out_params, ks0_raw);
}

/**
 * TLweKeySwitchKey destructor
 */
template<typename TORUS>
void TLweKeySwitchKey<TORUS>::destroy_obj(TLweKeySwitchKey<TORUS>* obj) {
    const int n = obj->n;
    const int t = obj->t;
    const int base = obj->base;
    del_obj_array(n*t*base, obj->ks0_raw);
    obj->~TLweKeySwitchKey();
}

