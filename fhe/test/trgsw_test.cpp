#include <include/gtest/gtest.h>
#include "../TRLwe.h"
#include "../TRGSW.h"

using namespace std;

TEST(TRGSW_TEST, trgsw_external_product) {


    int64_t N = 128;
    int64_t nblimbs_in = 2;


    int64_t alpha_bits = 80;
    int64_t bits_a = 15;
    UINT64 out_alpha_bits = alpha_bits - bits_a;

    BigTorusParams bt_params_in(nblimbs_in);
    TRGSWParams trgswParams(N, bt_params_in);

    shared_ptr<TLweKey> key = tlwe_keygen(trgswParams);

    TRGSW a(trgswParams);

    TRLwe b(trgswParams);

    TRLwe reps(trgswParams);


    int64_t *plaintext_a = new int64_t[a.params.N];

    for (UINT64 i = 0; i < a.params.N; i++) {
        plaintext_a[i] = (rand() % (1l << bits_a)) - (1l << (bits_a - 1));
    }

    intPoly_encrypt(a, plaintext_a, *key, alpha_bits);

    BigTorusPolynomial plaintext_b(N, bt_params_in);
    random(plaintext_b, nblimbs_in);

    native_encrypt(b, plaintext_b, *key, alpha_bits);

    external_product(reps, a, b, out_alpha_bits);

    BigTorusPolynomial phase(N, bt_params_in);

    native_phase(phase, reps, *key, alpha_bits);

    BigTorusPolynomial phase2(N, bt_params_in);
    fft_external_product(phase2, plaintext_a, plaintext_b, bits_a, limb_precision(out_alpha_bits));

    for (int64_t i = 0; i < N; i++) {
        ASSERT_LE(log2Diff(phase.getAT(i), phase2.getAT(i)), -1000);
    }
}

/*
TEST(TRGSW_TEST, trgsw_external_product_trivial) {
    int64_t N = 64;
    int64_t nblimbs = 6;
    int64_t ell;

    BigTorusParams bt_params(nblimbs);
    TLweParams tlwe_params_in(N, bt_params);
    TRGSWParams trgswParams(N, bt_params);

    TRLwe b(trgswParams);
    TRGSW a(trgswParams);

    //make a trivial ciphertext in b
    zero(b.a[0]);
    random(b.a[1], nblimbs);

    int64_t plaintextA = (rand() & 0xffff) - 32768;

    //make a trivial ciphertext in a
    a.fft_nlimbs = nblimbs;
    a.bits_a = 16;
    a.ell = ell;

    for (int j=0; j<2; j++) {
        for (int i=0; i<ell; i++) {
            for
        }
    }



}
*/
