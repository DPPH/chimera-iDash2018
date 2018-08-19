#ifndef FHE_TLWE_H
#define FHE_TLWE_H

#include <memory>
#include "BigTorusVector.h"

class TLweParams {
public:
    const UINT64 N;  ///< key size
    const BigTorusParams fixp_params; ///< fixp_params;

    TLweParams(const UINT64 N, const BigTorusParams &fixp_params);
};


/** A TLwe is an array of BigTorus elements, with fix point parameters */
class TLweKey {
public:
    int64_t *const key;

    NO_COPY(TLweKey);

    TLweKey(const TLweParams &params);

    ~TLweKey();
};


/** A TLwe is an array of BigTorus elements, with fix point parameters */
class TLwe : public BigTorusVector {
public:
    const TLweParams &params;

    TLwe(const TLweParams &params);

    ~TLwe();

    NO_COPY(TLwe);

    BigTorusRef getBT() const; ///< coef b as a BigTorus
    BigTorusRef getBT(); ///< coef b as a BigTorus
};

std::shared_ptr<TLweKey> tlwe_keygen(const TLweParams &params);

void zero(TLwe &tlwe);

void native_encrypt(TLwe &reps, const BigTorusRef &plaintext, const TLweKey &key, UINT64 alpha_bits = NA);

void native_phase(BigTorusRef reps, const TLwe &tlwe, const TLweKey &key, UINT64 alpha_bits = NA);

void slot_encrypt(TLwe &reps, const NTL::RR &plaintext, const TLweKey &key, UINT64 plaintext_precision);

NTL::RR slot_decrypt(const TLwe &tlwe, const TLweKey &key);

#endif //FHE_TLWE_H
