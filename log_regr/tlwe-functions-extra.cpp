#include "tlwe-functions-extra.h"

/**
 * fills the KeySwitching key array
 * @param result The (n x t x base) array of samples.
 *    result[i][j][k] encodes k.s[i]/base^(j+1)
 * @param out_key The TLwe key to encode all the output samples
 * @param out_alpha The standard deviation of all output samples
 * @param in_key The (binary) input Lwe key
 * @param n The size of the input key
 * @param t The precision of the keyswitch (technically, 1/2.base^t)
 * @param basebit Log_2 of base
 * @param N degree of the polynomials of output samples
 * @param k number of polynomials of output samples
 */
template<typename TORUS>
void TLweFunctionsExtra<TORUS>::CreateKeySwitchKey_fromArray(
    TLweSample<TORUS>** result,
    const TLweKey<TORUS>* out_key,
    const double out_alpha,
    const int* in_key,
    const int n,
    const int t,
    const int basebit)
{
    for(int i=0;i<n;i++) {
        for(int j=0;j<t;j++){
            TORUS x=in_key[i]*(TORUS(1)<<(TorusUtils<TORUS>::bit_cnt-(j+1)*basebit));
            TLweFunctions<TORUS>::SymEncryptT(&result[i][j],x,out_alpha,out_key);
            //printf("i,j,k,ki,x,phase=%d,%d,%d,%d,%d,%d\n",i,j,k,in_key->key[i],x,lwePhase(&result->ks[i][j][k],out_key));
        }
    }
}

/**
 * translates the message of the result sample by -sum(a[i].s[i]) where s is the secret
 * embedded in ks.
 * @param result the TLWE sample to translate by -sum(ai.si).
 * @param ks The (n x t x base) key switching key
 *    ks[i][j][k] encodes k.s[i]/base^(j+1)
 * @param params The common TLwe parameters of ks and result
 * @param ai The input torus array
 * @param n The size of the input key
 * @param t The precision of the keyswitch (technically, 1/2.base^t)
 * @param basebit Log_2 of base
 */
template<typename TORUS>
void TLweFunctionsExtra<TORUS>::KeySwitchTranslate_fromArray(
    TLweSample<TORUS>* result,
    const TLweSample<TORUS>** ks,
    const TLweParams<TORUS>* params,
    const TORUS* ai,
    const int n,
    const int t,
    const int basebit)
{
    const int base=1<<basebit;     // base=2 in [CGGI16]
    const TORUS prec_offset=TORUS(1)<<(TorusUtils<TORUS>::bit_cnt-(1+basebit*t)); //precision
    const int mask=base-1;

    typedef typename TorusUtils<TORUS>::UTORUS UTORUS;

    for (int i=0;i<n;i++){
        const UTORUS aibar=ai[i]+prec_offset;
        for (int j=0;j<t;j++){
            const UTORUS aij=(aibar>>(TorusUtils<TORUS>::bit_cnt-(j+1)*basebit)) & mask;
            TLweFunctions<TORUS>::SubMulTo(result,aij,&ks[i][j],params);
        }
    }
}

template<typename TORUS>
void TLweFunctionsExtra<TORUS>::CreateKeySwitchKey(
    TLweKeySwitchKey<TORUS>* result,
    const LweKey<TORUS>* in_key,
    const TLweKey<TORUS>* out_key)
{
    const int n=result->n;
    const int basebit=result->basebit;
    const int t=result->t;

    //TODO check the parameters

    CreateKeySwitchKey_fromArray(
        result->ks,
        out_key, out_key->params->alpha_min,
        in_key->key, n, t, basebit);
}

// sample=(a',b')
template<typename TORUS>
void TLweFunctionsExtra<TORUS>::KeySwitch(
    TLweSample<TORUS>* result,
    const TLweKeySwitchKey<TORUS>* ks,
    const LweSample<TORUS>* sample)
{
    const TLweParams<TORUS>* params=ks->out_params;
    const int n=ks->n;
    const int basebit=ks->basebit;
    const int t=ks->t;

    TLweFunctions<TORUS>::Clear(result, params);
    result->b->coefsT[0] = sample->b;

    TLweFunctionsExtra<TORUS>::KeySwitchTranslate_fromArray(
        result,
        (const TLweSample<TORUS>**) ks->ks, params,
        sample->a, n, t, basebit);
}
