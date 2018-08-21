#ifndef FHE_TLWE_H
#define FHE_TLWE_H

#include <memory>
#include "BigTorusVector.h"

#include "commons.h"
#include <cstdint>


class TLweParams {
public:
    const UINT64 N;  ///< key size
    const BigTorusParams fixp_params; ///< fixp_params;

    TLweParams(const UINT64 N, const BigTorusParams &fixp_params);
};

/** serialize:
 *  N:       UINT64
 *  param:   BigTorusParams
 */
void serializeTLweParams(std::ostream &out, const TLweParams &tLweParams);

/** serialize:
 *  N:       UINT64
 *  param:   BigTorusParams
 */
std::shared_ptr<TLweParams> deserializeTLweParams(std::istream &in);


/** A TLweKey  */
class TLweKey {
public:
    int64_t *const key;

    NO_COPY(TLweKey);

    TLweKey(const TLweParams &params);

    TLweKey(const UINT64 N);

    ~TLweKey();
};


/** serialize:
 *  magic number:   int64 on 8 bytes
 *  Key:     N*int8_t
 */
void serializeTLweKey(std::ostream &out, const TLweKey &key, const UINT64 N);

/** serialize:
 *  magic number:   int64 on 8 bytes
 *  key:     N*int8_t
 */
std::shared_ptr<TLweKey> deserializeTLweKey(std::istream &in, const UINT64 N);




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


/** serialize:

 *  value:        BigTorusVector
 */
void serializeTLwe(std::ostream &out, const TLwe &value);

/** serialize:
 *  content:        BigTorusVector
 */
std::shared_ptr<TLwe> deserializeTLwe(std::istream &in);




std::shared_ptr<TLweKey> tlwe_keygen(const TLweParams &params);

void zero(TLwe &tlwe);

void native_encrypt(TLwe &reps, const BigTorusRef &plaintext, const TLweKey &key, UINT64 alpha_bits = NA);

void native_phase(BigTorusRef reps, const TLwe &tlwe, const TLweKey &key, UINT64 alpha_bits = NA);

void slot_encrypt(TLwe &reps, const NTL::RR &plaintext, const TLweKey &key, UINT64 alpha_bits = NA);

NTL::RR slot_decrypt(const TLwe &tlwe, const TLweKey &key, UINT64 alpha_bits = NA);

#endif //FHE_TLWE_H
