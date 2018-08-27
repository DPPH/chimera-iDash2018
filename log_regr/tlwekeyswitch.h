#ifndef TLWEKEYSWITCH_H
#define TLWEKEYSWITCH_H

#include "tfhe_core.h"
#include "tfhe_io.h"

template<typename TORUS>
struct TLweKeySwitchKey: InitializerTag {
    int32_t n; ///< length of the input key: s'
    int32_t t; ///< decomposition length
    int32_t basebit; ///< log_2(base)
    int32_t base; ///< decomposition base: a power of 2
    const TLweParams<TORUS>* out_params; ///< params of the output key s
    TLweSample<TORUS>* ks0_raw;  //tableau qui contient tout les TLwe samples de taille n.t
    TLweSample<TORUS>** ks; ///< the keyswitch elements: a n.t matrix
    // de taille n pointe vers ks0 un tableau dont les cases sont espaceés de t positions

    TLweKeySwitchKey(int32_t n, int32_t t, int32_t basebit, const TLweParams<TORUS>* out_params, TLweSample<TORUS>* ks0_raw);
    ~TLweKeySwitchKey();
    TLweKeySwitchKey(const TLweKeySwitchKey<TORUS>&) = delete;
    void operator=(const TLweKeySwitchKey<TORUS>&) = delete;

    static void init_obj(
        TLweKeySwitchKey<TORUS>* obj,
        int32_t n,
        int32_t t,
        int32_t basebit,
        const TLweParams<TORUS>* out_params);
    static void destroy_obj(TLweKeySwitchKey<TORUS>* obj);

    void write(Ostream& out_stream) const;
    static TLweKeySwitchKey* read_new(Istream& inp_stream, const TLweParams<TORUS>* params);
};

// template struct TLweKeySwitchKey<Torus32>;
template struct TLweKeySwitchKey<Torus64>;

#endif //TLWEKEYSWITCH_H

