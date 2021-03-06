#include <include/gtest/gtest.h>
#include <NTL/RR.h>
#include "../TRLwe.h"
#include "../TRGSW.h"
#include "../BigFFT.h"
#include "../arithmetic.h"
#include "../BigTorus.h"

using namespace std;
using namespace NTL;

int64_t *debug_plaintext;
TLweKey *debug_key;

TEST(TRGSW_TEST, trgsw_external_product) {


    int64_t N = 4096;
    //int64_t n = 2 * N;
    int64_t nblimbs_in = 3;


    int64_t alpha_bits = 120; //signed
    int64_t bits_a = 15;
    UINT64 out_alpha_bits = alpha_bits - (32 + int(log2(N)));

    BigTorusParams bt_params_in(nblimbs_in);
    TRGSWParams trgswParams(N, bt_params_in);

    shared_ptr<TLweKey> key = tlwe_keygen(trgswParams);
    debug_key = key.get();
    //for (int64_t i = 0; i < N; i++) {
    //    key->key[i] = 0;
    //}

    TRGSW a(trgswParams);

    TRLwe b(trgswParams);

    TRLwe reps(trgswParams);


    int64_t *plaintext_a = new int64_t[a.params.N];
    debug_plaintext = plaintext_a;

    for (UINT64 i = 0; i < a.params.N; i++) {
        //plaintext_a[i] = ((i == 0) ? 1 : 0);
        plaintext_a[i] = (rand() % (1l << bits_a)) - (1l << (bits_a - 1));
    }

    intPoly_encrypt(a, plaintext_a, *key, alpha_bits);


    //verify the encryption of a
    BigTorusPolynomial phaseAC(N, bt_params_in);
    BigTorusPolynomial phaseRef(N, bt_params_in);

    //const BigComplex *powombar = fftAutoPrecomp.omegabar(n, nblimbs_in);

    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < int(a.ell); i++) {
            // the phase of a[j][i]
            native_phase_FFT(phaseAC, a.a[j][i][0], a.a[j][i][1], *key);
            // plaintext * phase(h_j,i)
            zero(phaseRef);
            phaseRef.getAT(0).limbs_end[-i / 2 - 1] = ((i % 2 == 0) ? (1ul << 32ul) : 1ul);
            if (j == 0) {
                mpn_neg(phaseRef.getAT(0).limbs_end - nblimbs_in,
                        phaseRef.getAT(0).limbs_end - nblimbs_in, nblimbs_in);
                fft_external_product(phaseRef, key->key, phaseRef, 1, nblimbs_in);
            }
            fft_external_product(phaseRef, plaintext_a, phaseRef, 16, nblimbs_in);
            //cout << "phase of a[" << j << "][" << i << "]: " << endl;
            for (int k = 0; k < N; k++) {
                //cout << j << "," << i << "," << k << ": "
                //     << log2Diff(phaseRef.getAT(k), phaseAC.getAT(k)) << endl;
                //cout << to_RR(phaseRef.getAT(k)) << endl;
                //cout << to_RR(phaseAC.getAT(k)) << endl;
                ASSERT_LE(log2Diff(phaseRef.getAT(k), phaseAC.getAT(k)), -alpha_bits + 2);
            }
        }
    }

    BigTorusPolynomial plaintext_b(N, bt_params_in);
    random(plaintext_b, nblimbs_in);

    native_encrypt(b, plaintext_b, *key, alpha_bits);

    cout << "external 1 start: " << clock() / double(CLOCKS_PER_SEC) << endl;
    external_product(reps, a, b, out_alpha_bits);
    cout << "external 2 start: " << clock() / double(CLOCKS_PER_SEC) << endl;
    external_product(reps, a, b, out_alpha_bits);
    cout << "external 2 end: " << clock() / double(CLOCKS_PER_SEC) << endl;

    BigTorusPolynomial phase(N, bt_params_in);

    native_phase(phase, reps, *key, alpha_bits);

    BigTorusPolynomial phase2(N, bt_params_in);
    fft_external_product(phase2, plaintext_a, plaintext_b, bits_a, limb_precision(out_alpha_bits));

    for (int64_t i = 0; i < N; i++) {
        ASSERT_LE(log2Diff(phase.getAT(i), phase2.getAT(i)), -out_alpha_bits);
    }
}


TEST(TRGSW_BLINDROTATE_TEST, trgsw_blind_rotate) {
    //int64_t N = 4096;
    int64_t N = 256;
    //int64_t n_in = 500;
    int64_t n_in = 50;
    int64_t nblimbs = 2;
    int64_t alpha_bits = 120; //signed
    int64_t out_alpha_bits = alpha_bits - (32 + int(log2(N)));

    BigTorusParams bt_params(nblimbs);

    TRGSWParams trgswParams(N, bt_params);

    BigTorusPolynomial phase(N, bt_params);

    TRLwe reps(trgswParams);
    TRGSW *c = new_TRGSW_array(n_in, trgswParams);
    shared_ptr<TLweKey> key = tlwe_keygen(trgswParams);
    int64_t b = (rand() % (2 * N));


    int64_t *s = new int64_t[n_in];
    int64_t *a = new int64_t[n_in];
    int64_t power = -b;

    cout << "start encrypt: " << clock() / double(CLOCKS_PER_SEC) << endl;
    for (int i = 0; i < n_in; i++) {
        s[i] = rand() % 2;
        a[i] = rand() % (2 * N);
        power += s[i] * a[i];
    }
#pragma omp parallel for
    for (int i = 0; i < n_in; i++) {
        int_encrypt(c[i], s[i], *key, alpha_bits);
    }
    cout << "end encrypt: " << clock() / double(CLOCKS_PER_SEC) << endl;



    //make a trivial ciphertext in reps
    zero(reps.a[0]);
    random(reps.a[1], nblimbs);
    TRLwe phase_ref(trgswParams);
    copy(phase_ref, reps);


    cout << "start blind rotate: " << clock() / double(CLOCKS_PER_SEC) << endl;
    blind_rotate(reps, b, a, c, n_in, out_alpha_bits);
    cout << "end blind rotate: " << clock() / double(CLOCKS_PER_SEC) << endl;

    native_phase(phase, reps, *key, alpha_bits);
    rotate(phase_ref, phase_ref, power);


    for (int64_t i = 0; i < N; i++) {
        EXPECT_LE(log2Diff(phase.getAT(i), phase_ref.a[1].getAT(i)), -out_alpha_bits + 2);
    }
    delete[] a;
    delete[] s;
    delete_TRGSW_array(n_in, c);
}

TEST(TRGSW_TEST, trlwe_internal_product) {

    //int64_t N = 4096;
    int64_t N = 128;

    int64_t L_a = 110; //90; //level expo of a
    int64_t L_b = 80; //90; //level expo of b
    int64_t rho = 18; //precision bits
    int64_t Ns2 = N / 2;
    int64_t n = N * 2;
    int64_t nblimbs_a = limb_precision(L_a + rho + 5);
    int64_t nblimbs_b = limb_precision(L_b + rho + 5);
    int64_t L = std::min(L_a, L_b) - rho - 1; //output level expo
    int64_t nblimbs_reps = limb_precision(L + rho + 5);
    int64_t alpha_rk = L + rho + 32 + log2(N / 2);
    int64_t nblimbs_rk = limb_precision(alpha_rk);

    BigTorusParams bt_params_a(nblimbs_a, 0, L_a);
    BigTorusParams bt_params_b(nblimbs_b, 0, L_b);
    BigTorusParams bt_params_reps(nblimbs_reps, 0, L);

    TRLweParams trlweParams_a(N, bt_params_a);
    TRLweParams trlweParams_b(N, bt_params_b);
    TRLweParams trlweParams_reps(N, bt_params_reps);

    BigTorusParams bt_params_rk(nblimbs_rk);
    TRGSWParams trgswParams_rk(N, bt_params_rk);

    shared_ptr<TLweKey> key = tlwe_keygen(trlweParams_reps);
    debug_key = key.get();

    TRLwe a(trlweParams_a);
    TRLwe b(trlweParams_b);
    TRLwe reps(trlweParams_reps);
    TRGSW rk(trgswParams_rk);

    //here, we draw complex slots at random (w.r.t. their respective level) and compute the product
    BigComplex *slots_a = new_BigComplex_array(Ns2, nblimbs_a);
    BigComplex *slots_b = new_BigComplex_array(Ns2, nblimbs_b);
    BigComplex *slots_prod = new_BigComplex_array(Ns2, nblimbs_reps);
    BigReal *coefs_a = new_BigReal_array(N, nblimbs_a);
    BigReal *coefs_b = new_BigReal_array(N, nblimbs_b);
    BigReal *coefs_prod = new_BigReal_array(N, nblimbs_reps);

    RR::SetPrecision(nblimbs_a * BITS_PER_LIMBS + 32);
    for (int i = 0; i < Ns2; i++) {
        RR are = (random_RR() - to_RR(0.5));
        RR aim = (random_RR() - to_RR(0.5));
        RR bre = (random_RR() - to_RR(0.5));
        RR bim = (random_RR() - to_RR(0.5));
        to_BigReal(slots_a[i].real, are * power2_RR(-L_a));
        to_BigReal(slots_a[i].imag, aim * power2_RR(-L_a));
        to_BigReal(slots_b[i].real, bre * power2_RR(-L_b));
        to_BigReal(slots_b[i].imag, bim * power2_RR(-L_b));
        RR rre = are * bre - aim * bim;
        RR rim = aim * bre + are * bim;
        to_BigReal(slots_prod[i].real, rre * power2_RR(-L));
        to_BigReal(slots_prod[i].imag, rim * power2_RR(-L));
    }
    //Then, we compute their FFT and consider a and b as plaintext
    FFT(coefs_a, slots_a, n, fftAutoPrecomp.omegabar(n, nblimbs_a));
    FFT(coefs_b, slots_b, n, fftAutoPrecomp.omegabar(n, nblimbs_b));
    FFT(coefs_prod, slots_prod, n, fftAutoPrecomp.omegabar(n, nblimbs_reps));

    BigTorusPolynomial plaintext_a(N, bt_params_a);
    BigTorusPolynomial plaintext_b(N, bt_params_b);
    for (int i = 0; i < N; i++) {
        shift_toBigTorus(plaintext_a.getAT(i), coefs_a[i], 0);
        shift_toBigTorus(plaintext_b.getAT(i), coefs_b[i], 0);
    }

    //we encrypt a, b, and the relin key, and perform a heaan product
    native_encrypt(a, plaintext_a, *key, L_a + rho + 1);
    native_encrypt(b, plaintext_b, *key, L_b + rho + 1);

    intPoly_encrypt(rk, key->key, *key, alpha_rk);

    cout << "time: " << clock() / double(CLOCKS_PER_SEC) << endl;
    fixp_internal_product(reps, a, b, rk, rho);
    cout << "time: " << clock() / double(CLOCKS_PER_SEC) << endl;
    fixp_internal_product(reps, a, b, rk, rho);
    cout << "time: " << clock() / double(CLOCKS_PER_SEC) << endl;

    // finally, we verify that the product matches the expected product
    BigTorusPolynomial phase(N, bt_params_reps);
    native_phase(phase, reps, *key, L + rho);

    for (int64_t i = 0; i < N; i++) {
        EXPECT_LE(log2Diff(to_RR(phase.getAT(i)), to_RR(coefs_prod[i])), -L - rho + 5);
    }

    delete_BigReal_array(N, coefs_a);
    delete_BigReal_array(N, coefs_b);
    delete_BigReal_array(N, coefs_prod);
    delete_BigComplex_array(Ns2, slots_a);
    delete_BigComplex_array(Ns2, slots_b);
    delete_BigComplex_array(Ns2, slots_prod);
}

TEST(TRGSW_TEST, trlwe_internal_product_square) {

    //int64_t N = 4096;
    int64_t N = 512;

    int64_t L_a = 21; //90; //level expo of a
    int64_t L_b = 21; //90; //level expo of b
    int64_t rho = 14; //precision bits
    int64_t Ns2 = N / 2;
    int64_t n = N * 2;
    int64_t nblimbs_a = limb_precision(L_a + rho + 5);
    int64_t nblimbs_b = limb_precision(L_b + rho + 5);
    int64_t L = 1; //std::min(L_a, L_b) - rho; //output level expo
    int64_t nblimbs_reps = limb_precision(L + rho + 5);
    int64_t alpha_rk = L + rho + 32 + log2(N / 2);
    int64_t nblimbs_rk = limb_precision(alpha_rk);

    BigTorusParams bt_params_a(nblimbs_a, 0, L_a);
    BigTorusParams bt_params_b(nblimbs_b, 0, L_b);
    BigTorusParams bt_params_reps(nblimbs_reps, 0, L);

    TRLweParams trlweParams_a(N, bt_params_a);
    TRLweParams trlweParams_b(N, bt_params_b);
    TRLweParams trlweParams_reps(N, bt_params_reps);

    BigTorusParams bt_params_rk(nblimbs_rk);
    TRGSWParams trgswParams_rk(N, bt_params_rk);

    shared_ptr<TLweKey> key = tlwe_keygen(trlweParams_reps);
    debug_key = key.get();

    TRLwe a(trlweParams_a);
    TRLwe b(trlweParams_b);
    TRLwe reps(trlweParams_reps);
    TRGSW rk(trgswParams_rk);

    //here, we draw complex slots at random (w.r.t. their respective level) and compute the product
    BigComplex *slots_a = new_BigComplex_array(Ns2, nblimbs_a);
    BigComplex *slots_b = new_BigComplex_array(Ns2, nblimbs_b);
    BigComplex *slots_prod = new_BigComplex_array(Ns2, nblimbs_reps);
    BigComplex *slots_prod2 = new_BigComplex_array(Ns2, nblimbs_reps);
    BigReal *coefs_a = new_BigReal_array(N, nblimbs_a);
    BigReal *coefs_b = new_BigReal_array(N, nblimbs_b);
    BigReal *coefs_prod = new_BigReal_array(N, nblimbs_reps); //TODO
    BigReal *coefs_prod2 = new_BigReal_array(N, nblimbs_reps); //TODO

    RR::SetPrecision(nblimbs_a * BITS_PER_LIMBS + 32);
    for (int i = 0; i < Ns2; i++) {
        RR are = (random_RR() - to_RR(0.5));
        RR aim = (random_RR() - to_RR(0.5));
        RR bre = are; // (random_RR() - to_RR(0.5));
        RR bim = aim; // (random_RR() - to_RR(0.5));
        to_BigReal(slots_a[i].real, are * power2_RR(-L_a));
        to_BigReal(slots_a[i].imag, aim * power2_RR(-L_a));
        to_BigReal(slots_b[i].real, bre * power2_RR(-L_b));
        to_BigReal(slots_b[i].imag, bim * power2_RR(-L_b));
        RR rre = are * bre - aim * bim;
        RR rim = aim * bre + are * bim;
        to_BigReal(slots_prod[i].real, rre * power2_RR(-L));
        to_BigReal(slots_prod[i].imag, rim * power2_RR(-L));
    }
    //Then, we compute their FFT and consider a and b as plaintext
    FFT(coefs_a, slots_a, n, fftAutoPrecomp.omegabar(n, nblimbs_a));
    FFT(coefs_b, slots_b, n, fftAutoPrecomp.omegabar(n, nblimbs_b));

    BigTorusPolynomial plaintext_a(N, bt_params_a);
    BigTorusPolynomial plaintext_b(N, bt_params_b);
    for (int i = 0; i < N; i++) {
        shift_toBigTorus(plaintext_a.getAT(i), coefs_a[i], 0);
        shift_toBigTorus(plaintext_b.getAT(i), coefs_b[i], 0);
    }

    //we encrypt a, b, and the relin key, and perform a heaan product
    native_encrypt(a, plaintext_a, *key, L_a + rho + 1);
    native_encrypt(b, plaintext_b, *key, L_b + rho + 1);

    intPoly_encrypt(rk, key->key, *key, alpha_rk);

    cout << "time: " << clock() / double(CLOCKS_PER_SEC) << endl;
    fixp_internal_product(reps, a, b, rk, rho);
    cout << "time: " << clock() / double(CLOCKS_PER_SEC) << endl;
    fixp_internal_product(reps, a, b, rk, rho);
    cout << "time: " << clock() / double(CLOCKS_PER_SEC) << endl;

    // finally, we verify that the product matches the expected product
    BigTorusPolynomial phase(N, bt_params_reps);
    native_phase(phase, reps, *key, L + rho);
    iFFT(slots_prod2, phase);


    FFT(coefs_prod, slots_prod, n, fftAutoPrecomp.omegabar(n, nblimbs_reps)); //TODO
    FFT(coefs_prod2, slots_prod2, n, fftAutoPrecomp.omegabar(n, nblimbs_reps)); //coefs from phase

    RR diffRe;
    diffRe = 0;
    for (int64_t i = 0; i < N; i++) {
        EXPECT_LE(log2Diff(to_RR(phase.getAT(i)), to_RR(coefs_prod[i])), -L - rho + log2(sqrt(N)));
        EXPECT_LE(log2Diff(to_RR(coefs_prod[i]), to_RR(coefs_prod2[i])), -L - rho + log2(sqrt(N)));
        diffRe += sqr(to_RR(phase.getAT(i)) - to_RR(coefs_prod[i]));
    }
    RR diffCpl;
    diffCpl = 0;
    for (int64_t i = 0; i < N / 2; i++) {
        if (i == 0) {
            cout << "re: " << to_RR(slots_prod[i].real) << endl;
            cout << "im: " << to_RR(slots_prod[i].imag) << endl;
            cout << "re2: " << to_RR(slots_prod2[i].real) << endl;
            cout << "im2: " << to_RR(slots_prod2[i].imag) << endl;
        }
        EXPECT_LE(log2Diff(to_RR(slots_prod[i].real), to_RR(slots_prod2[i].real)), -L - rho + 5 + log2(sqrt(N)));
        EXPECT_LE(log2Diff(to_RR(slots_prod[i].imag), to_RR(slots_prod2[i].imag)), -L - rho + 5 + log2(sqrt(N)));
        diffCpl += sqr(to_RR(slots_prod[i].real) - to_RR(slots_prod2[i].real))
                   + sqr(to_RR(slots_prod[i].imag) - to_RR(slots_prod2[i].imag));
    }
    cout << "DiffRe: " << diffRe << endl;
    cout << "DiffCpl: " << diffCpl << endl;

    delete_BigReal_array(N, coefs_a);
    delete_BigReal_array(N, coefs_b);
    delete_BigReal_array(N, coefs_prod);
    delete_BigReal_array(N, coefs_prod2);
    delete_BigComplex_array(Ns2, slots_a);
    delete_BigComplex_array(Ns2, slots_b);
    delete_BigComplex_array(Ns2, slots_prod);
    delete_BigComplex_array(Ns2, slots_prod2);
}

TEST(TRGSW_TEST, trlwe_internal_product_const) {

    int64_t N = 4096;
    //int64_t N = 128;

    int64_t L_a = 100; //level expo of a
    int64_t L_b = 110; //level expo of b
    int64_t rho = 20; //precision bits
    int64_t nblimbs_a = limb_precision(L_a + rho + 5);
    int64_t nblimbs_b = limb_precision(L_b + rho + 5);
    int64_t L = std::min(L_a, L_b) - rho - 1; //output level expo
    int64_t nblimbs_reps = limb_precision(L + rho + 5);
    int64_t alpha_rk = L + rho + 32 + log2(N / 2);
    int64_t nblimbs_rk = limb_precision(alpha_rk);

    BigTorusParams bt_params_a(nblimbs_a, 0, L_a);
    BigTorusParams bt_params_b(nblimbs_b, 0, L_b);
    BigTorusParams bt_params_reps(nblimbs_reps, 0, L);

    TRLweParams trlweParams_a(N, bt_params_a);
    TRLweParams trlweParams_b(N, bt_params_b);
    TRLweParams trlweParams_reps(N, bt_params_reps);

    BigTorusParams bt_params_rk(nblimbs_rk);
    TRGSWParams trgswParams_rk(N, bt_params_rk);

    shared_ptr<TLweKey> key = tlwe_keygen(trlweParams_reps);
    debug_key = key.get();

    TRLwe a(trlweParams_a);
    TRLwe b(trlweParams_b);
    TRLwe reps(trlweParams_reps);
    TRGSW rk(trgswParams_rk);

    BigTorusPolynomial plaintext_a(N, bt_params_a);
    BigTorusPolynomial plaintext_b(N, bt_params_b);
    int64_t p_a = (rand() % (1ul << rho)) - (1ul << (rho - 1));
    int64_t p_b = (rand() % (1ul << rho)) - (1ul << (rho - 1));

    zero(plaintext_a);
    zero(plaintext_b);

    RR::SetPrecision(alpha_rk + 10);
    RR ra;
    RR rb;
    ra = to_RR(long(p_a)) / power2_RR(L_a + rho);
    rb = to_RR(long(p_b)) / power2_RR(L_b + rho);

    cout << "plain_a" << p_a / power2_RR(rho) << endl;
    cout << "plain_b" << p_b / power2_RR(rho) << endl;
    cout << "plain_prod" << p_a * p_b / power2_RR(2 * rho) << endl;

    to_torus(plaintext_a.getAT(0), ra);
    to_torus(plaintext_b.getAT(0), rb);

    native_encrypt(a, plaintext_a, *key, L_a + rho + 1);
    native_encrypt(b, plaintext_b, *key, L_b + rho + 1);

    intPoly_encrypt(rk, key->key, *key, alpha_rk);

    cout << "time: " << clock() / double(CLOCKS_PER_SEC) << endl;
    fixp_internal_product(reps, a, b, rk, rho);
    cout << "time: " << clock() / double(CLOCKS_PER_SEC) << endl;
    fixp_internal_product(reps, a, b, rk, rho);
    cout << "time: " << clock() / double(CLOCKS_PER_SEC) << endl;

    BigTorusPolynomial phase(N, bt_params_reps);

    native_phase(phase, reps, *key, L + rho);


    RR exp_phase = to_RR(long(p_a * p_b)) / power2_RR(L + 2 * rho);

    for (int64_t i = 1; i < N; i++) {
        EXPECT_LE(log2Diff(to_RR(phase.getAT(i)), to_RR(0)), -L - rho + 5);
    }
    EXPECT_LE(log2Diff(to_RR(phase.getAT(0)), exp_phase), -L - rho + 5);
}
