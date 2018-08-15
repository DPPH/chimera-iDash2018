#include <include/gtest/gtest.h>
#include "../TRLwe.h"
#include "../TRGSW.h"

using namespace std;

TEST(TRGSW_TEST, trgsw_external_product) {


    int64_t N = 128;
    int64_t nblimbs_in = 2;


    int64_t alpha_bits = 80;
    UINT64 out_alpha_bits = alpha_bits - 15;


    BigTorusParams bt_params_in(nblimbs_in);

    TLweParams tlwe_params_in(N, bt_params_in);
    shared_ptr<TLweKey> key = tlwe_keygen(tlwe_params_in);

    TRGSWParams trgswParams(N, bt_params_in);
    TRGSW a(trgswParams);


    TRLweParams trlweParams(N, bt_params_in);
    TRLwe b(trlweParams);

    TRLwe reps(trlweParams);


    int64_t *plaintext_a = new int64_t[a.params.N];

    for (UINT64 i = 0; i < a.params.N; i++) {
        plaintext_a[i] = (rand() % 65536) - 32768;
    }

    intPoly_encrypt(a, plaintext_a, *key, alpha_bits);

    BigTorusPolynomial plaintext_b(N, bt_params_in);
    random(plaintext_b, nblimbs_in);

    native_encrypt(b, plaintext_b, *key, alpha_bits);

    external_product(reps, a, b, out_alpha_bits);

    BigTorusPolynomial phase(N, bt_params_in);

    native_phase(phase, reps, *key, alpha_bits);

    BigTorusPolynomial phase2(N, bt_params_in);
    fft_external_product(phase2, plaintext_a, plaintext_b, 15, limb_precision(out_alpha_bits));

    for (UINT64 i = 0; i < N; i++) {
        ASSERT_LE(log2Diff(phase.getAT(i), phase2.getAT(i)), -1000);
    }
}
