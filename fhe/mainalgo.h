#ifndef FHE_MAINALGO_H
#define FHE_MAINALGO_H

#include "TRGSW.h"
#include <NTL/mat_RR.h>

class TRGSWMatrix {
public:
    TRGSW** data;
    int rows;
    int cols;
private:
    TRGSWMatrix()
};

/**
 * encrypt S as packed heaan TRGSW slots
 *
 */
std::shared_ptr<TRGSW> encrypt_S(
        NTL::mat_RR plaintext,
        int64_t N, int64_t alpha_bits, int64_t plaintext_precision_bits) {

}

#endif //FHE_MAINALGO_H
