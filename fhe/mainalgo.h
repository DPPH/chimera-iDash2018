#ifndef FHE_MAINALGO_H
#define FHE_MAINALGO_H

#include "TRGSW.h"
#include "TRLwe.h"
#include <NTL/mat_RR.h>

static const int default_plaintext_precision = 18;

class TRGSWMatrix {
public:
    TRGSW **const data;
    const int64_t rows;
    const int64_t cols;
private:
    TRGSW *const data_raw;
public:
    TRGSWMatrix(int64_t rows, int64_t cols, const TRGSWParams &params) :
            data(new TRGSW *[rows]),
            rows(rows),
            cols(cols),
            data_raw(new_TRGSW_array(rows * cols, params)) {
        for (int i = 0; i < rows; ++i) data[i] = data_raw + i * cols;
    }

    ~TRGSWMatrix() {
        delete_TRGSW_array(rows * cols, data_raw);
        delete[] data;
    }

    NO_COPY(TRGSWMatrix);
};

class TRLweMatrix {
public:
    TRLwe **const data;
    const int64_t rows;
    const int64_t cols;
private:
    TRLwe *const data_raw;
public:
    TRLweMatrix(int64_t rows, int64_t cols, const TRLweParams &params) :
            data(new TRLwe *[rows]),
            rows(rows),
            cols(cols),
            data_raw(new_TRLwe_array(rows * cols, params)) {
        for (int i = 0; i < rows; ++i) data[i] = data_raw + i * cols;
    }

    ~TRLweMatrix() {
        delete_TRLwe_array(rows * cols, data_raw);
        delete[] data;
    }

    NO_COPY(TRLweMatrix);
};

class TRLWEVector {
public:
    TRLwe *const data;
    const int64_t length;
public:
    TRLWEVector(int64_t length, const TRLweParams &params) :
            data(new_TRLwe_array(length, params)),
            length(length) {
    }

    ~TRLWEVector() {
        delete_TRLwe_array(length, data);
    }

    NO_COPY(TRLWEVector);
};

/**
 * encrypt S as packed heaan TRGSW slots
 *
 */
std::shared_ptr<TRGSWMatrix>
encrypt_S(NTL::mat_RR plaintext, const TLweKey &key, int64_t N, int64_t alpha_bits,
          int64_t plaintext_precision_bits = default_plaintext_precision);

/**
 * encrypt X as ind heaan TRGSW
 *
 */
std::shared_ptr<TRGSWMatrix>
encrypt_X(NTL::mat_RR plaintext, const TLweKey &key, int64_t N, int64_t alpha_bits,
          int64_t plaintext_precision_bits = default_plaintext_precision);

/**
 * encrypt y and/or p as individual TRLWE slots
 */
std::shared_ptr<TRLWEVector>
encrypt_individual_trlwe(NTL::vec_RR plaintext, const TLweKey &key, int64_t N,
                         int64_t level_expo,
                         int64_t plaintext_expo = NA,
                         int64_t plaintext_precision_bits = default_plaintext_precision);

/**
  * @brief encrypted vector - Matrix multiplication -- autodeduce parameters
  *
  * @param v  the vector to multiply, encrypted as individual TRLWE
  * @param A  the matrix to multiply, encrypted as packed TRGSW slots
  * @param target_level_expo        if available, the targeted output level exponent
  * @param override_plaintext_expo  if available, the plaintext exponent of the result (overrides the default)
  * @param plaintext_precision_bits the plaintext precision of the operation
  * @return
  */
std::shared_ptr<TRLWEVector>
mat_vec_prod(const TRLWEVector &v, const TRGSWMatrix &A,
             int64_t target_level_expo = NA,
             int64_t override_plaintext_exponent = NA,
             int64_t plaintext_precision_bits = default_plaintext_precision);

/**
 * @brief decrypt TRLWE slots as a vec<RR>
 */
NTL::vec_RR decrypt_heaan_packed_trlwe(const TRLWEVector &ciphertext, const TLweKey &key, int64_t length);

/**
 * @brief decrypt TRLWE slots as a vec<RR>
 */
NTL::vec_RR decrypt_individual_trlwe(const TRLWEVector &ciphertext, const TLweKey &key, int64_t length);


/**
 * @brief multiply two TRLWE vectors
 */
std::shared_ptr<TRLWEVector>
product_ind_TRLWE(const TRLWEVector &a, const TRLWEVector &b, const TRGSW &rk, int64_t target_level_expo = NA,
                  int64_t override_plaintext_exponent = NA,
                  int64_t plaintext_precision_bits = default_plaintext_precision);


/**
 * @brief substract two TRLWE vectors
 */
std::shared_ptr<TRLWEVector>
substract_ind_TRLWE(const TRLWEVector &a, const TRLWEVector &b, int64_t target_level_expo = NA,
                    int64_t override_plaintext_exponent = NA,
                    int64_t plaintext_precision_bits = default_plaintext_precision);


/**
 * @brief add two TRLWE vectors
 */
std::shared_ptr<TRLWEVector>
add_ind_TRLWE(const TRLWEVector &a, const TRLWEVector &b, int64_t target_level_expo = NA,
              int64_t override_plaintext_exponent = NA,
              int64_t plaintext_precision_bits = default_plaintext_precision);
/**
 * @brief copy a to the output
 */
std::shared_ptr<TRLWEVector> copy(const TRLWEVector &a);

/**
 * @brief copy -a to the output
 */
std::shared_ptr<TRLWEVector> neg(const TRLWEVector &a);


/**
 * @brief Compute w= p(1-p)
 */
std::shared_ptr<TRLWEVector> compute_w(const TRLWEVector &p, const TRGSW &rk, int64_t target_level_expo = NA,
                                       int64_t override_plaintext_exponent = NA,
                                       int64_t plaintext_precision_bits = default_plaintext_precision);

/**
 * @brief Compute A= X^t* S * W
 * A  matrix of TRLWE k+1 to l
 * X  matrix of TRGSW n to k+1
 * S  matrix of TRGSW n to l
 * W vector of n TRLWE
 *
 */
std::shared_ptr<TRLweMatrix>
compute_A(const TRGSWMatrix &X, const TRGSWMatrix &S, const TRLWEVector &W,
          int64_t target_level_expo = NA,
          int64_t override_plaintext_exponent = NA,
          int64_t plaintext_precision_bits = default_plaintext_precision);





#endif //FHE_MAINALGO_H
