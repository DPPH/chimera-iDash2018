#ifndef FHE_TLWE_H
#define FHE_TLWE_H

#include <memory>
#include "BigFixP.h"

class TLweParams {
public:
    const uint64_t N;  ///< key size
    const BigFixPParams fixp_params; ///< fixp_params;

    TLweParams(const uint64_t N, const BigFixPParams &fixp_params);
};


/** A TLwe is an array of BigTorus elements, with fix point parameters */
class TLweKey {
public:
    int8_t *const key;
    const TLweParams &params;

    NO_COPY(TLweKey);

    TLweKey(const TLweParams &params);

    ~TLweKey();
};


/** A TLwe is an array of BigTorus elements, with fix point parameters */
class TLwe {
public:
    uint64_t *const limbs;
    const TLweParams &params;

    NO_COPY(TLwe);

    TLwe(const TLweParams &params);

    ~TLwe();

    BigTorusRef getAT(uint64_t i) const; ///< coef a[i] as a BigTorus
    BigTorusRef getAT(uint64_t i); ///< coef a[i] as a BigTorus
    BigTorusRef getBT() const; ///< coef b as a BigTorus
    BigTorusRef getBT(); ///< coef b as a BigTorus
    BigFixPRef getAF(uint64_t i) const; ///< coef a[i] as a BigFixedPoint
    BigFixPRef getAF(uint64_t i); ///< coef a[i] as a BigFixedPoint
    BigFixPRef getBF() const; ///< coef b as a BigFixedPoint
    BigFixPRef getBF(); ///< coef b as a BigFixedPoint
};

std::shared_ptr<TLweKey> tlwe_keygen(const TLweParams &params);

void zero(TLwe &tlwe);

void native_encrypt(TLwe &reps, const BigTorusRef &plaintext, const TLweKey &key, uint64_t alpha_bits = NA);

void native_phase(BigTorusRef reps, const TLwe &tlwe, const TLweKey &key, uint64_t alpha_bits = NA);

void slot_encrypt(TLwe &reps, const NTL::RR &plaintext, const TLweKey &key, uint64_t alpha_bits = NA);

NTL::RR slot_decrypt(const TLwe &tlwe, const TLweKey &key, uint64_t alpha_bits = NA);

#endif //FHE_TLWE_H
