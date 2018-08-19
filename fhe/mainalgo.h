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
 * encrypt y and/or p as individual TRLWE slots
 */
std::shared_ptr<TRLWEVector>
encrypt_individual_trlwe(NTL::vec_RR plaintext, const TLweKey &key, int64_t N,
                         int64_t level_expo,
                         int64_t plaintext_expo = NA,
                         int64_t plaintext_precision_bits = default_plaintext_precision);

/**
 * @brief encrypted vector - Matrix multiplication -- autodeduce parameters
 */
std::shared_ptr<TRLWEVector>
mat_vec_prod(const TRLWEVector &v, const TRGSWMatrix &A,
             int64_t target_level_expo = NA, int64_t plaintext_precision_bits = default_plaintext_precision);

/**
 * @brief decrypt TRLWE slots as RR
 */
NTL::vec_RR decrypt_slots(const TRLWEVector &ciphertext, const TLweKey &key, int64_t length);


#endif //FHE_MAINALGO_H
