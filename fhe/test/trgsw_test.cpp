#include <include/gtest/gtest.h>
#include <NTL/RR.h>
#include "../TRLwe.h"
#include "../TRGSW.h"
#include "../BigFFT.h"
#include "../arithmetic.h"

using namespace std;
using namespace NTL;

TEST(TRGSW_TEST, trgsw_external_product) {


    int64_t N = 128;
    int64_t n = 2 * N;
    int64_t nblimbs_in = 3;


    int64_t alpha_bits = 112;
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
    //verify the encryption of a
    BigReal *aa = new_BigReal_array(N, nblimbs_in);
    BigReal *ab = new_BigReal_array(N, nblimbs_in);
    TRLwe ac(trgswParams);
    BigTorusPolynomial phaseAC(N, bt_params_in);

    const BigComplex *powombar = fftAutoPrecomp.omegabar(n, nblimbs_in);

    for (int i = 0; i < int(a.ell); i++) {
        //the phase of a[1][i] should be (plaintext, 0)/2^32.(i+1)
        FFT(aa, a.a[1][i][0], n, powombar);
        FFT(ab, a.a[1][i][1], n, powombar);
        for (int k = 0; k < N; k++) {
            to_BigTorus(ac.a[0].getAT(k), aa[k], 0, nblimbs_in);
        }
        for (int k = 0; k < N; k++) {
            to_BigTorus(ac.a[1].getAT(k), ab[k], 0, nblimbs_in);
        }
        //compute the phase of ac
        native_phase(phaseAC, ac, *key, alpha_bits);
        //compare with plaintext / 2^32.(i+1)
        for (int k = 0; k < N; k++) {
            RR::SetPrecision(nblimbs_in * BITS_PER_LIMBS);
            RR actual_k = to_RR(phaseAC.getAT(k));
            RR expected_k = to_RR(plaintext_a[k]) / pow(to_RR(2), to_RR(32 * (i + 1)));
            //cout << "-----" << endl;
            //cout << actual_k << endl;
            //cout << expected_k << endl;
            ASSERT_LE(log2Diff(actual_k, expected_k), -alpha_bits + 1);
        }
    }

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

    delete_BigReal_array(N, aa);
    delete_BigReal_array(N, ab);
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
