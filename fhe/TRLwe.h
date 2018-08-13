#ifndef FHE_TRLWE_H
#define FHE_TRLWE_H


#include "BigTorusPolynomial.h"
#include "TLwe.h"

class TRLweParams : public TLweParams {
public:
    TRLweParams(const UINT64 N, const BigTorusParams &fixp_params);

    ~TRLweParams() = default;
};

class TRLwe {
public:
    const TRLweParams &params;
    BigTorusPolynomial *const a;

    TRLwe(const TRLweParams &params);

    ~TRLwe();
};

/** @brief allocates a matrix of TRLWE */
TRLwe **new_TRLweMatrix(UINT64 rows, UINT64 cols, const TRLweParams &params);

/** @brief frees a matrix of TRLWE allocated with new_TRLWEMatrix */
void delete_TRLweMatrix(TRLwe **data, UINT64 rows, UINT64 cols, const TRLweParams &params);


struct pubKsKey {
    static const UINT64 BgBits = 128;
    const TLweParams &in_params;
    const TRLweParams &out_params;
    const UINT64 l_dec;

    TRLwe **const kskey;
    BigTorus bitDecomp_in_offset; // sum Bg/2 Bg^i
    __int128 bitDecomp_out_offset; // -Bg/2

    /**
     * constructs a pubKsKey structure, ksKey and offsets are not initialized
     * @param in_params
     * @param out_params
     * @param l_dec decomposition length in base 2^128
     */
    pubKsKey(const TLweParams &in_params, const TRLweParams &out_params, const UINT64 l_dec);

    ~pubKsKey();
};



/**
 * @brief generation of PublicKeySwitch key
 */
std::shared_ptr<pubKsKey> ks_keygen(const TRLweParams &out_params, const TLweParams &in_params, const TLweKey &in_key, const TLweKey &out_key, const UINT64 out_limb_prec);


void pubKS(TRLwe &out, TLwe &in, pubKsKey &ks, const UINT64 out_prec_limbs);


/**
 * @brief create a trivial ciphertext of a constant
 */
void trivial(TRLwe &out, const BigTorusRef &in, const UINT64 out_limb_prec);

/** @brief get the j-th coef of the base 2^128 decomposition between [0 and 2^128[
 * the output can be signed, we don't care */
__int128 bitdecomp_coef128(const BigTorusRef &tmpDec, UINT64 j, const UINT64 limb_prec);


/** @brief out = out - aij * in
 *  WARNING: input limb precision = out_limb_prec + 2 */
void subMul(TRLwe &out, __int128 aij, const TRLwe &in, const UINT64 out_limb_prec);

void native_encrypt(TRLwe &reps, const BigTorusPolynomial &plaintext, const TLweKey &key, UINT64 alpha_bits);

void native_phase(BigTorusPolynomial &reps, const TRLwe &c, const TLweKey &key, UINT64 alpha_bits);

void zero(TRLwe &out);

#endif //FHE_TRLWE_H
