#ifndef FHE_MAINALGO_H
#define FHE_MAINALGO_H

#include "TRGSW.h"
#include <NTL/mat_RR.h>

class TRGSWMatrix {
public:
    TRGSW** data;
    int64_t rows;
    int64_t cols;
private:
    TRGSW *data_raw;
public:
    TRGSWMatrix(int64_t rows, int64_t cols, const TRGSWParams &params);

    ~TRGSWMatrix();

    NO_COPY(TRGSWMatrix);
};

/**
 * encrypt S as packed heaan TRGSW slots
 *
 */
std::shared_ptr<TRGSWMatrix>
encrypt_S(NTL::mat_RR plaintext, const TLweKey &key, int64_t N, int64_t alpha_bits, int64_t plaintext_precision_bits);


#endif //FHE_MAINALGO_H
