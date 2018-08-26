#ifndef FHE_TRGSW_H
#define FHE_TRGSW_H


#include "BigTorusPolynomial.h"
#include "TRLwe.h"

class TRGSWParams : public TRLweParams {
public:
    static const UINT64 Bgbits = 32;
    static const UINT64 max_ell = 5;

    TRGSWParams(const UINT64 N, const BigTorusParams &fixp_params);

    ~TRGSWParams() = default;
};

/** serialize:
 *  param:   TRGSWParams
 */
void serializeTRGSWParams(std::ostream &out, const TRGSWParams &params);

/** serialize:
 *  param:   TRGSWParams
 */
std::shared_ptr<TRGSWParams> deserializeTRGSWParams(std::istream &in);


class TRGSW {
public:
    const TRGSWParams &params;
    BigComplex *(a[2][TRGSWParams::max_ell][2]);
    UINT64 plaintext_exponent; //the exponent we should multiply the trgsw to obtain the original plaintext
    UINT64 bits_a;             //bits of the norm 1 of the plaintext
    UINT64 fft_nlimbs; //fft limbs in all BigComplex
    UINT64 ell;        //actual decomposition length

    TRGSW(const TRGSWParams &params);

    ~TRGSW();
};

/** serialize:
 *  magic number:   int64 on 8 bytes
 *  a:        tab[][][] of BigComplex
 *  plaintext_exponent: UINT64
 *  bits_a: UINT64
 *  fft_nlimbs: UINT64
 *  ell: UINT64
 */
void serializeTRGSWContent(std::ostream &out, const TRGSW &value);

/** serialize:
 *  magic number:   int64 on 8 bytes
 *  a:        tab[][][] of BigComplex
 *  plaintext_exponent: UINT64
 *  bits_a: UINT64
 *  fft_nlimbs: UINT64
 *  ell: UINT64
 */
void deserializeTRGSWContent(std::istream &in, TRGSW &reps);


/** serialize:
 * params:
 * value:
 */
void serializeTRGSW(std::ostream &out, const TRGSW &value);

/** serialize:
 * params:
 * value:
 */
std::shared_ptr<TRGSW> deserializeTRGSW(std::istream &in);


void intPoly_encrypt(TRGSW &reps, const int64_t *plaintext, const TLweKey &key, UINT64 alpha_bits);

void int_encrypt(TRGSW &reps, const int64_t plaintext, const TLweKey &key, UINT64 alpha_bits);

void external_product(TRLwe &reps, const TRGSW &a, const TRLwe &b, int64_t out_alpha_bits);

void
native_phase_FFT(BigTorusPolynomial &reps, BigComplex *a, BigComplex *b, const TLweKey &key, const int64_t lshift = 0);

// reps = b + c.(a-b)
void cmux(TRLwe &reps, const TRGSW &c, const TRLwe &a, const TRLwe &b, int64_t out_alpha_bits);

// reps *= X^power
void rotate(TRLwe &out, const TRLwe &in, int64_t power);

void rotate_diff(TRLwe &out, const TRLwe &in, int64_t power);

//reps *= X^-(b-sum c_i a_i)
void blind_rotate(TRLwe &reps, int64_t b, int64_t *a, const TRGSW *c, int64_t n_in, int64_t out_alpha_bits);


/** create an array of TRGSW */
TRGSW *new_TRGSW_array(UINT64 size, const TRGSWParams &params);

/** delete an array of TRGSW */
void delete_TRGSW_array(UINT64 size, TRGSW *array);


void fixp_internal_product(TRLwe &reps, const TRLwe &a, const TRLwe &b, const TRGSW &rk, int precision_bits);


#endif //FHE_TRGSW_H
