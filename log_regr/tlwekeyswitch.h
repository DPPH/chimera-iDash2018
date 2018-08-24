#define TLWEKEYSWITCH_H

#include "tfhe_core.h"
#include "tfhe_alloc.h"
#include "lweparams.h"
#include "lwesamples.h"
#include "tlwe.h"

template<typename TORUS>
    struct TLweKeySwitchKey: InitializerTag {
       int n; ///< length of the input key: s'
       int t; ///< decomposition length
       int basebit; ///< log_2(base)
       int base; ///< decomposition base: a power of 2
       const TLweParams<TORUS>* out_params; ///< params of the output key s
       TLweSample<TORUS>* ks0_raw;
       //tableau qui contient tout les TLwe samples de taille nlbase
       TLweSample<TORUS>** ks1_raw; // de taille nl  pointe vers un tableau ks0_raw dont les cases sont espaceés de base positions
       TLweSample<TORUS>*** ks;
       ///< the keyswitch elements: a n.l.base matrix
       // de taille n pointe vers ks1 un tableau dont les cases sont espaceés de ell positions

        TLweKeySwitchKey(int n, int t, int basebit, const TLweParams<TORUS>* out_params, TLweSample<TORUS>* ks0_raw);
        virtual ~TLweKeySwitchKey() = 0;
        TLweKeySwitchKey(const TLweKeySwitchKey<TORUS>&) = delete;
        void operator=(const TLweKeySwitchKey<TORUS>&) = delete;

        static void init_obj(
                TLweKeySwitchKey<TORUS>* obj,
                 int n,
                 int t,
                int basebit,
                const TLweParams<TORUS>* out_params);
        static void destroy_obj(TLweKeySwitchKey<TORUS>* obj);
    };


template struct TLweKeySwitchKey<Torus32>;
template struct TLweKeySwitchKey<Torus64>;

//allocates and initialize the LweKeySwitchKey structure
//(equivalent of the C++ new)
template<typename TORUS>
inline TLweKeySwitchKey<TORUS>* new_TLweKeySwitchKey(int n, int t, int basebit, const TLweParams<TORUS>* out_params) {
    return new_obj<LweKeySwitchKey<TORUS>>(n, t, basebit, out_params);
}

//might not use it
template<typename TORUS>
inline TLweKeySwitchKey<TORUS>* new_TLweKeySwitchKey_array(int nbelts, int n, int t, int basebit, const TLweParams<TORUS>* out_params) {
    return new_obj_array<TLweKeySwitchKey<TORUS>>(nbelts, n, t, basebit, out_params);
}

//destroys and frees the TLweKeySwitchKey structure
//(equivalent of the C++ delete)
template<typename TORUS>
inline void delete_TLweKeySwitchKey(TLweKeySwitchKey<TORUS>* obj) {
    del_obj<TLweKeySwitchKey<TORUS>>(obj);
}

#endif //TLWEKEYSWITCH_H

