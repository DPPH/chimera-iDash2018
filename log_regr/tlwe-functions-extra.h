#ifndef TLWE_FUNCTIONS_EXTRA_H
#define TLWE_FUNCTIONS_EXTRA_H

///@file
///@brief This file contains the operations on Lwe samples

#include "tlwekeyswitch.h"

template<typename TORUS>
struct TLweFunctionsExtra
{
private:
    static void CreateKeySwitchKey_fromArray(
        TLweSample<TORUS>** result,
        const TLweKey<TORUS>* out_key,
        const double out_alpha,
        const int* in_key,
        const int n,
        const int t,
        const int basebit);

    static void KeySwitchTranslate_fromArray(
        TLweSample<TORUS>* result,
        const TLweSample<TORUS>** ks,
        const TLweParams<TORUS>* params,
        const TORUS* ai,
        const int n,
        const int t,
        const int basebit);

public:
    /**
     * creates a Key Switching Key between the two keys
     */
    static void CreateKeySwitchKey(
        TLweKeySwitchKey<TORUS>* result,
        const LweKey<TORUS>* in_key,
        const TLweKey<TORUS>* out_key);

    /**
     * applies keySwitching
     */
    static void KeySwitch(
        TLweSample<TORUS>* result,
        const TLweKeySwitchKey<TORUS>* ks,
        const LweSample<TORUS>* sample);
};

// template struct TLweFunctionsExtra<Torus32>;
template struct TLweFunctionsExtra<Torus64>;

#endif //TLWE_FUNCTIONS_EXTRA_H
