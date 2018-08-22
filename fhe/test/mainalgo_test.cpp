#include <include/gtest/gtest.h>
#include "../mainalgo.h"
#include "../arithmetic.h"
#include <NTL/mat_RR.h>

using namespace std;
using namespace NTL;

TEST(MAINALGO, encrypt_S) {
    int rows = 1; //should be 300 in reality
    int full_cols = 1000; //should be 10000 in reality
    int N = 256;  // 4096 in reality
    int plaintext_precision = 18;
    int alpha_bits = 128;
    int Ns2 = N / 2;
    int packed_cols = (full_cols + Ns2 - 1) / Ns2;

    mat_RR plaintext(INIT_SIZE, rows, full_cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < full_cols; j++) {
            plaintext[i][j] = rand() % 2;
        }
    }
    BigTorusParams bigTorusParams(0);
    TLweParams blo(N, bigTorusParams); //TODO keep only N to initialize a key
    shared_ptr<TLweKey> key = tlwe_keygen(blo);
    shared_ptr<TRGSWMatrix> c = encrypt_S(plaintext, *key, N, alpha_bits, plaintext_precision);

    BigComplex *slot_test = new_BigComplex_array(Ns2, c->data[0][0].fft_nlimbs);
    BigTorusParams g_params(c->data[0][0].fft_nlimbs);
    BigTorusPolynomial poly_test(N, g_params);

    ASSERT_EQ(c->rows, rows);
    ASSERT_EQ(c->cols, packed_cols);
    //line [1][0] of trgsw should encrypt plaintext*(2^bits_a-pe) / 2^32
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < packed_cols; j++) {
            native_phase_FFT(poly_test, c->data[i][j].a[1][0][0], c->data[i][j].a[1][0][1], *key);
            iFFT(slot_test, poly_test);
            for (int k = 0; k < Ns2; k++) {
                if (j * Ns2 + k >= plaintext.NumCols()) break;
                //slot_test[k].real should be equal to plaintext[k]*(2^bits_a-pe-32)
                //slot_test[k].imag should be equal to 0
                //cerr << to_RR(poly_test.getAT(k)) << endl;
                //cerr << to_RR(slot_test[k].imag) << endl;
                //cerr << to_RR(slot_test[k].real) << endl;
                //cerr << plaintext[i][j*Ns2+k]*power2_RR(plaintext_precision-c->data[0][0].plaintext_exponent-32) << endl;
                RR actual = to_RR(slot_test[k].real) * power2_RR(-plaintext_precision + 32);
                RR expext = plaintext[i][j * Ns2 + k] * power2_RR(-c->data[0][0].plaintext_exponent);
                EXPECT_LE(log2Diff(to_RR(slot_test[k].imag), to_RR(0)), -64);
                EXPECT_LE(log2Diff(actual, expext), -plaintext_precision + log2(N) / 2 + 1);
            }
        }
    }
    delete_BigComplex_array(Ns2, slot_test);
}


TEST(MAINALGO, product_ind_TRLWE) {
    int64_t N = 128;
    int length = 10;

    int64_t L_a = 80; //level expo of a
    int64_t L_b = 110; //level expo of b
    int64_t rho = 20; //precision bits


    vec_RR vec_a;
    vec_a.SetLength(length);
    vec_RR vec_b;
    vec_b.SetLength(length);

    for (int i = 0; i < length; i++) {
        vec_a[i] = random_RR();
        vec_b[i] = random_RR();
    }


    BigTorusParams bt_params_key(0, 0, 0);
    TRLweParams trlweParams_key(N, bt_params_key);
    shared_ptr<TLweKey> key = tlwe_keygen(trlweParams_key);


    shared_ptr<TRLWEVector> a = encrypt_individual_trlwe(vec_a, *key, N, L_a);
    shared_ptr<TRLWEVector> b = encrypt_individual_trlwe(vec_b, *key, N, L_b);


    int64_t alpha_rk = 128;
    int64_t nblimbs_rk = limb_precision(alpha_rk);
    BigTorusParams bt_params_rk(nblimbs_rk);
    TRGSWParams trgswParams_rk(N, bt_params_rk);
    TRGSW rk(trgswParams_rk);
    intPoly_encrypt(rk, key->key, *key, alpha_rk);


    shared_ptr<TRLWEVector> resp = product_ind_TRLWE(*a, *b, rk, 70);

    vec_RR target_resp = decrypt_individual_trlwe(*resp, *key, length);

    for (int i = 0; i < length; i++) {
        EXPECT_LE(log2Diff(target_resp[i], vec_a[i] * vec_b[i]), -rho);

    }


}


TEST(MAINALGO, substract_ind_TRLWE) {

}
