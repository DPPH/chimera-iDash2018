#include <include/gtest/gtest.h>
#include "../mainalgo.h"
#include <NTL/mat_RR.h>

using namespace std;
using namespace NTL;

TEST(MAINALGO, encrypt_S) {
    int rows = 31;
    int cols = 300;
    int N = 128;
    int plaintext_precision = 16;
    int alpha_bits = 128;

    mat_RR plaintext(INIT_SIZE, rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            plaintext[i][j] = rand() % 2;
        }
    }
    BigTorusParams bigTorusParams(0);
    TLweParams blo(N, bigTorusParams); //TODO keep only N to initialize a key
    shared_ptr<TLweKey> key = tlwe_keygen(blo);
    shared_ptr<TRGSWMatrix> c = encrypt_S(plaintext, *key, N, alpha_bits, plaintext_precision);

    ASSERT_EQ(c->rows, rows);
    ASSERT_EQ(c->cols, ceil(cols / double(N / 2)));


}


TEST(MAINALGO, product_ind_TRLWE) {

}


TEST(MAINALGO, substract_ind_TRLWE) {

}
