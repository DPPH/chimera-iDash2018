#include <include/gtest/gtest.h>
#include <NTL/mat_RR.h>
#include "../BigReal.h"
#include "../BigComplex.h"
#include "../BigFFT.h"
#include "../BigTorusPolynomial.h"

NTL_CLIENT

double log2Diff(const RR &a, const RR &b);

TEST(FFT_TEST, precompute_iFFT) {
    int nblimbs = 2;
    int n = 1024;
    //int N = n / 2;
    //int ns4 = n / 4;

    BigComplex *powomega = precomp_iFFT(n, nblimbs);


    RR::SetPrecision(nblimbs * BITS_PER_LIMBS);
    RR _2pi = 2 * ComputePi_RR();
    for (int i = 0; i < n; i++) {
        RR angle = to_RR(i) * _2pi / to_RR(n);
        RR tcos = cos(angle);
        RR xcos = to_RR(powomega[i].real);
        //cerr << i << " tcos:" << tcos <<  " xcos:" << xcos << endl;
        ASSERT_LE(log2Diff(tcos, xcos), -nblimbs * BITS_PER_LIMBS + 3);
        RR tsin = sin(angle);
        RR xsin = to_RR(powomega[i].imag);
        ASSERT_LE(log2Diff(tsin, xsin), -nblimbs * BITS_PER_LIMBS + 3);
    }

    clear_precomp_iFFT(powomega);
}

TEST(FFT_TEST, precompute_FFT) {
    int nblimbs = 2;
    int n = 1024;
    //int N = n / 2;
    //int ns4 = n / 4;

    BigComplex *powombar = precomp_FFT(n, nblimbs);


    RR::SetPrecision(nblimbs * BITS_PER_LIMBS);
    RR _2pi = 2 * ComputePi_RR();
    for (int i = 0; i < n; i++) {
        RR angle = to_RR(-i) * _2pi / to_RR(n);
        RR tcos = cos(angle);
        RR xcos = to_RR(powombar[i].real);
        //cerr << i << " tcos:" << tcos <<  " xcos:" << xcos << endl;
        ASSERT_LE(log2Diff(tcos, xcos), -nblimbs * BITS_PER_LIMBS + 3);
        RR tsin = sin(angle);
        RR xsin = to_RR(powombar[i].imag);
        ASSERT_LE(log2Diff(tsin, xsin), -nblimbs * BITS_PER_LIMBS + 3);
    }
    clear_precomp_iFFT(powombar);
}

TEST(FFT_TEST, bijection_iFFT_FFT) {
    int nblimbs = 3;
    int n = 512;
    int N = n / 2;
    int ns4 = n / 4;
    int64_t PREC = nblimbs * BITS_PER_LIMBS;

    BigComplex *powombar = precomp_FFT(n, nblimbs);
    BigComplex *powomega = precomp_iFFT(n, nblimbs);

    //create an array of bigreal input
    BigReal *in = new_BigReal_array(N, nblimbs);
    //same for the output
    BigComplex *out = new_BigComplex_array(ns4, nblimbs);
    //create an array of bigreal input
    BigReal *out2 = new_BigReal_array(N, nblimbs);

    RR::SetPrecision(PREC);
    vec_RR test_in;
    test_in.SetLength(N);
    // fill random values
    for (int i = 0; i < N; i++) {
        test_in[i] = random_RR() / ns4 / ns4;
        to_BigReal(in[i], test_in[i]);
    }
    // iFFT
    iFFT(out, in, n, powomega);
    // FFT
    FFT(out2, out, n, powombar);
    //test equality
    for (int i = 0; i < N; i++) {
        RR out2_i = to_RR(out2[i]);
        ASSERT_LE(log2Diff(out2_i, test_in[i]), -PREC + 3);
    }

    //clean up
    delete_BigReal_array(N, out2);
    delete_BigComplex_array(ns4, out);
    delete_BigReal_array(N, in);
    clear_precomp_FFT(powombar);
    clear_precomp_iFFT(powomega);
}

TEST(FFT_TEST, bijection_FFT_iFFT) {
    int nblimbs = 2;
    int n = 256;
    int N = n / 2;
    int ns4 = n / 4;
    int64_t PREC = nblimbs * BITS_PER_LIMBS;

    BigComplex *powombar = precomp_FFT(n, nblimbs);
    BigComplex *powomega = precomp_iFFT(n, nblimbs);

    //create an array of bigreal input
    BigComplex *in = new_BigComplex_array(ns4, nblimbs);
    //same for the output
    BigReal *out = new_BigReal_array(N, nblimbs);
    //create an array of bigreal input
    BigComplex *out2 = new_BigComplex_array(ns4, nblimbs);

    RR::SetPrecision(PREC);
    vec_RR test_in;
    test_in.SetLength(N);
    // fill random values
    for (int i = 0; i < ns4; i++) {
        test_in[2 * i] = random_RR() / ns4 / ns4;
        test_in[2 * i + 1] = random_RR() / ns4 / ns4;
        to_BigReal(in[i].real, test_in[2 * i]);
        to_BigReal(in[i].imag, test_in[2 * i + 1]);
    }
    // FFT
    FFT(out, in, n, powombar);
    // iFFT
    iFFT(out2, out, n, powomega);
    //test equality
    for (int i = 0; i < ns4; i++) {
        RR out2_i_r = to_RR(out2[i].real);
        RR out2_i_i = to_RR(out2[i].imag);
        //cerr << i << endl;
        //here, the precision is less good because of the rshift in FFT
        ASSERT_LE(log2Diff(out2_i_r, test_in[2 * i]), -PREC + 20);
        ASSERT_LE(log2Diff(out2_i_i, test_in[2 * i + 1]), -PREC + 20);
    }

    //clean up
    delete_BigComplex_array(ns4, out2);
    delete_BigReal_array(N, out);
    delete_BigComplex_array(ns4, in);

    clear_precomp_FFT(powombar);
    clear_precomp_iFFT(powomega);
}

TEST(FFT_EXT_MUL_TEST, external_product_FFT) {
    for (int64_t N : {32, 64, 256}) {
        for (int64_t nblimbs_in : {1, 2, 3}) {
            int64_t nblimbs_out = nblimbs_in;

            BigTorusParams params_in(nblimbs_in);
            BigTorusParams params_out(nblimbs_out);

            BigTorusPolynomial b(N, params_in);

            BigTorusPolynomial out(N, params_out);
            BigTorusPolynomial out1(N, params_out);

            random(b, nblimbs_in);

            int64_t bits_a = 32;
            int64_t _2am1m1 = (1 << (bits_a - 1)) - 1;
            int64_t _2am1 = (1 << (bits_a)) - 1;

            int64_t *a = new int64_t[N];

            for (int64_t i = 0; i < N; i++) {
                a[i] = (rand() & _2am1) - _2am1m1;
            }
            //cout << clock() / double(CLOCKS_PER_SEC) << " avant fft " << endl;

            fft_external_product(out, a, b, bits_a, nblimbs_out);

            //cout << clock() / double(CLOCKS_PER_SEC) << " avant fft " << endl;

            fft_external_product(out, a, b, bits_a, nblimbs_out);

            //cout << clock() / double(CLOCKS_PER_SEC) << " apres fft " << endl;

            naive_external_product(out1, a, b, nblimbs_out);

            //cout << clock() / double(CLOCKS_PER_SEC) << " apres ext product " << endl;
            for (int64_t i = 0; i < N; i++) {
                //cout << "i: " << i << endl;
                //RR::SetPrecision(nblimbs_out * BITS_PER_LIMBS);
                //RR::SetOutputPrecision(nblimbs_out * BITS_PER_LIMBS / log2(10.));
                //cout << "out_i : " << to_RR(out.getAT(i)) << endl;
                //cout << "out1_i: " << to_RR(out1.getAT(i)) << endl;
                ASSERT_LE(log2Diff(out.getAT(i), out1.getAT(i)), -nblimbs_out * BITS_PER_LIMBS);
            }

            delete[] a;
        }
    }
}


TEST(FFT_INTERNAL_MUL_TEST, internal_product_FFT) {

    int64_t N = 8;
    int64_t nblimbs_in = 5;
    int64_t nblimbs_out = 5;

    BigTorusParams params_in(nblimbs_in);
    BigTorusParams params_out(nblimbs_out);

    BigTorusPolynomial b(N, params_in);
    BigTorusPolynomial a(N, params_in);

    BigTorusPolynomial out(N, params_out);
    BigTorusPolynomial out1(N, params_out);

    random(b, nblimbs_in);
    random(a, nblimbs_in);

    fft_internal_product(out, a, b, nblimbs_out);
    //TODO avec quoi comparer
    //naive_external_product(out1, a, b, nblimbs_out);
    // for (UINT64 i = 0; i < N; i++) {
    //    ASSERT_LE(log2Diff(out.getAT(i), out1.getAT(i)), -nblimbs_out * BITS_PER_LIMBS);

    //}
}

TEST(FFT_BIGREALPOLY_PRODUCT, fft_BigRealPoly_product) {

    for (int64_t N : {16, 64, 128}) {
        for (int64_t fft_nlimbs : {1, 3, 5, 6}) {

            BigReal *a = new_BigReal_array(N, fft_nlimbs);
            BigReal *b = new_BigReal_array(N, fft_nlimbs);
            RR::SetPrecision(fft_nlimbs * BITS_PER_LIMBS);


            RR *ra = new RR[N];
            RR *rb = new RR[N];
            RR *rreps = new RR[N];
            for (int i = 0; i < N; i++) {
                ra[i] = random_RR();
                to_BigReal(a[i], ra[i]);
                rb[i] = random_RR();
                to_BigReal(b[i], rb[i]);
            }

            BigReal *reps = new_BigReal_array(N, fft_nlimbs);
            fft_BigRealPoly_product(reps, a, b, N, fft_nlimbs);

            for (int64_t i = 0; i < N; i++) {
                RR ri;
                ri = 0;
                for (int64_t j = 0; j <= i; j++) {
                    ri += ra[j] * rb[i - j];

                }
                for (int64_t j = i + 1; j < N; j++) {
                    ri -= ra[j] * rb[N + i - j];
                }
                rreps[i] = ri;
            }


            for (int i = 0; i < N; ++i) {
                EXPECT_LE(log2Diff(to_RR(reps[i]), rreps[i]), -fft_nlimbs * BITS_PER_LIMBS + log2(N) + 2);
            }
            delete_BigReal_array(N, a);
            delete_BigReal_array(N, b);
            delete_BigReal_array(N, reps);
            delete[] ra;
            delete[] rb;
            delete[] rreps;
        }
    }
}

