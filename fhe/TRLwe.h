#ifndef FHE_TRLWE_H
#define FHE_TRLWE_H


#include "BigTorusPolynomial.h"
#include "TLwe.h"
#include <memory>
#include "commons.h"
#include <cstdint>
#include <NTL/mat_RR.h>


class TRLweParams : public TLweParams {
public:
    TRLweParams(const UINT64 N, const BigTorusParams &fixp_params);

    ~TRLweParams() = default;
};


/** serialize:
 *  magic number:   int64 on 8 bytes
 *  value:        TRLweParams
 */
void serializeTRLweParams(std::ostream &out, const TRLweParams &value);

/** serialize:
 *  magic number:   int64 on 8 bytes
 *  content:        TRLweParams
 */
std::shared_ptr<TRLweParams> deserializeTRLweParams(std::istream &in);




class TRLwe {
public:
    const TRLweParams &params;
    BigTorusPolynomial *const a;

    TRLwe(const TRLweParams &params);

    ~TRLwe();
};


/** serialize:
 *  magic number:   int64 on 8 bytes
 *  value:        N* size(TLWE)
 */
void serializeTRLweContent(std::ostream &out, const TRLwe &value);

/** serialize:
 *  magic number:   int64 on 8 bytes
 *  value:        N* size(TLWE)
 */
void deserializeTRLweContent(std::istream &in, TRLwe &reps);



/** serialize:
 *  magic number:   int64 on 8 bytes
 *  value:        N* size(TLWE)
 */
void serializeTRLwe(std::ostream &out, const TRLwe &value);

/** serialize:
 *   magic number:   int64 on 8 bytes
 *  content:        N* size(TLWE)
 */
std::shared_ptr<TRLwe> deserializeTRLwe(std::istream &in);


/** @brief allocates a matrix of TRLWE */
TRLwe **new_TRLweMatrix(UINT64 rows, UINT64 cols, const TRLweParams &params);

/** @brief frees a matrix of TRLWE allocated with new_TRLWEMatrix */
void delete_TRLweMatrix(TRLwe **data, UINT64 rows, UINT64 cols, const TRLweParams &params);


struct pubKsKey128 {
    static const UINT64 BgBits = 128;
    const TLweParams &in_params;
    const TRLweParams &out_params;
    const TRLweParams ks_params;
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
    pubKsKey128(const TLweParams &in_params, const TRLweParams &out_params, TRLweParams &&ks_params,
                const UINT64 l_dec);

    ~pubKsKey128();
};

struct pubKsKey32 {
    static const UINT64 BgBits = 32;
    const TLweParams &in_params;
    const TRLweParams &out_params;
    const TRLweParams ks_params;
    const UINT64 l_dec;

    TRLwe **const kskey;
    BigTorus bitDecomp_in_offset; // sum Bg/2 Bg^i
    int64_t bitDecomp_out_offset; // -Bg/2

    /**
     * constructs a pubKsKey structure, ksKey and offsets are not initialized
     * @param in_params
     * @param out_params
     * @param l_dec decomposition length in base 2^128
     */
    pubKsKey32(const TLweParams &in_params, const TRLweParams &out_params, TRLweParams &&ks_params, const UINT64 l_dec);

    ~pubKsKey32();
};

/** serialize:
 *  magic number:   int64 on 8 bytes
 *  BgBits:         UINT64
 *  in_params:     TLweParams
 *  out_params:    TRLweParams:
 *  ks_params:     TRLweParams
 *  l_dec:         UINT64
 *  kskey:         tab of TRLwe
 *  bitDecomp_in_offset: BigTorus
 *  bitDecomp_out_offset: int64_t
 */
void serializepubKsKey32(std::ostream &out, const pubKsKey32 &key);

/** serialize:
 *  magic number:   int64 on 8 bytes
 *  BgBits:         UINT64
 *  in_params:     TLweParams
 *  out_params:    TRLweParams:
 *  ks_params:     TRLweParams
 *  l_dec:         UINT64
 *  kskey:         tab of TRLwe
 *  bitDecomp_in_offset: BigTorus
 *  bitDecomp_out_offset: int64_t
 */
std::shared_ptr<pubKsKey32> deserializepubKsKey32(std::istream &in);

/**
 * @brief generation of PublicKeySwitch key
 * @param out_alpha_bits alpha of the output of the keyswitch (expect the output to be 128-bit more noisy)
 */
std::shared_ptr<pubKsKey128>
ks_keygen128(const TRLweParams &out_params, const TLweParams &in_params,
             const TLweKey &in_key, const TLweKey &out_key,
             const UINT64 out_alpha_bits);

/**
 * @brief generation of PublicKeySwitch key
 * @param out_alpha_bits alpha of the output of the keyswitch (expect the output to be 128-bit more noisy)
 */
std::shared_ptr<pubKsKey32>
ks_keygen32(const TRLweParams &out_params, const TLweParams &in_params,
            const TLweKey &in_key, const TLweKey &out_key,
            const UINT64 out_alpha_bits);


void pubKS128(TRLwe &out, const TLwe &in, const pubKsKey128 &ks, const UINT64 out_prec_limbs);

void pubKS32(TRLwe &out, const TLwe &in, const pubKsKey32 &ks, const UINT64 out_prec_limbs);


/**
 * @brief create a trivial ciphertext of a constant
 */
void trivial(TRLwe &out, const BigTorusRef &in, const UINT64 out_limb_prec);

/** @brief get the j-th coef of the base 2^128 decomposition between [0 and 2^128[
 * the output can be signed, we don't care */
__int128 bitdecomp_coef128(const BigTorusRef &tmpDec, UINT64 j, const UINT64 limb_prec);

/** @brief add the signed bitdecomp in offset.
 * This function must be called once on a Torus before doing a bit decomposition */
void bitdecomp_signed_offset32_apply(BigTorusRef reps, const BigTorusRef &source);

/** @brief get the j-th coef of the base 2^128 decomposition between [0 and 2^32[
 * the output is signed */
int64_t bitdecomp_coef32(const BigTorusRef &tmpDec, UINT64 position, const UINT64 limb_prec);

/** @brief get the j-th coef of the base 2^128 decomposition between [-2^31 and 2^31[
 * the output is signed */
int64_t bitdecomp_signed_coef32(const BigTorusRef &tmpDec, UINT64 position, const UINT64 limb_prec);



/** @brief out = out - aij * in
 *  WARNING: input limb precision = out_limb_prec + 2 */
void subMul128(TRLwe &out, __int128 aij, const TRLwe &in, const UINT64 out_limb_prec);

/** @brief out = out - aij * in
 *  WARNING: input limb precision = out_limb_prec + 2 */
void subMul64(TRLwe &out, int64_t aij, const TRLwe &in, const UINT64 out_limb_prec);

void native_encrypt(TRLwe &reps, const BigTorusPolynomial &plaintext, const TLweKey &key, UINT64 alpha_bits);

void native_phase(BigTorusPolynomial &reps, const TRLwe &c, const TLweKey &key, UINT64 alpha_bits);

void fixp_encrypt(TRLwe &reps, const NTL::vec_RR &plaintext, const TLweKey &key, UINT64 plaintext_precision);

void fixp_trivial(TRLwe &reps, const NTL::vec_RR &plaintext, UINT64 plaintext_precision);

NTL::vec_RR fixp_decrypt(const TRLwe &tlwe, const TLweKey &key);

/** @brief decrypt the N/2 complex slots inside the ciphertext (N/2 real first, N/2 imag in the end) */
NTL::vec_RR slot_decrypt(const TRLwe &tlwe, const TLweKey &key);

void fixp_encrypt_number(TRLwe &reps, const NTL::RR &plaintext, const TLweKey &key, UINT64 plaintext_precision);

void fixp_trivial_number(TRLwe &reps, const NTL::RR &plaintext, UINT64 plaintext_precision);

NTL::RR fixp_decrypt_number(const TRLwe &tlwe, const TLweKey &key);

void zero(TRLwe &out);

void add(TRLwe &reps, const TRLwe &a, const TRLwe &b);

void sub(TRLwe &reps, const TRLwe &a, const TRLwe &b);

void neg(TRLwe &out, const TRLwe &in);

void copy(TRLwe &out, const TRLwe &in);

/** @brief left shift the all torus by exactly shift_bits */
void lshift(TRLwe &out, const TRLwe &in, int64_t shift_bits);

/** create an array of TRGSW */
TRLwe *new_TRLwe_array(UINT64 size, const TRLweParams &params);

/** delete an array of TRGSW */
void delete_TRLwe_array(UINT64 size, TRLwe *array);

void fixp_add(TRLwe &reps, const TRLwe &a, const TRLwe &b, UINT64 out_precision_bits = NA);

void fixp_sub(TRLwe &reps, const TRLwe &a, const TRLwe &b, UINT64 out_precision_bits = NA);

void fixp_public_product(TRLwe &reps, const TRLwe &a, double value);

#endif //FHE_TRLWE_H

