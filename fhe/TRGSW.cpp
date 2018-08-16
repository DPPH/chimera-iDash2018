#include <cassert>
#include "TRGSW.h"
#include "BigFFT.h"

#include <NTL/RR.h>

NTL_CLIENT;

double log2Diff(const RR &a, const RR &b);


TRGSW::TRGSW(const TRGSWParams &params) : params(params), bits_a(0), fft_nlimbs(0), ell(0) {
    for (UINT64 k = 0; k < 2; k++) {
        for (UINT64 i = 0; i < params.max_ell; i++) {
            for (UINT64 j = 0; j < 2; j++) {
                a[k][i][j] = nullptr;
            }
        }
    }
}

TRGSW::~TRGSW() {
    for (UINT64 k = 0; k < 2; k++) {
        for (UINT64 i = 0; i < params.max_ell; i++) {
            for (UINT64 j = 0; j < 2; j++) {
                if (a[k][i][j])
                    delete_BigComplex_array(params.N / 2, a[k][i][j]);
            }
        }
    }
}


TRGSWParams::TRGSWParams(const UINT64 N, const BigTorusParams &fixp_params) :
        TRLweParams(N, fixp_params) {}


void binary_encrypt(TRGSW &reps, const UINT64 plaintext, const TLweKey &key, UINT64 alpha_bits) {

    int64_t *poly_plaintext = new int64_t[reps.params.N];
    poly_plaintext[0] = plaintext;
    for (UINT64 i = 1; i < reps.params.N; i++) {
        poly_plaintext[i] = 0;
    }
    intPoly_encrypt(reps, poly_plaintext, key, alpha_bits);
}

void intPoly_encrypt(
        TRGSW &reps, const int64_t *plaintext, const TLweKey &key, UINT64 alpha_bits) {
    const UINT64 N = reps.params.N;
    const UINT64 n = N * 2;
    const UINT64 Ns2 = N / 2;
    const UINT64 nb_bits_decomp = alpha_bits;
    const UINT64 ell = UINT64(ceil(double(nb_bits_decomp) / double(reps.params.Bgbits)));
    assert_dramatically(ell <= reps.params.max_ell, "TRGSW params ell is too small");
    const UINT64 fft_prec_bits = nb_bits_decomp + 32 + (UINT64) log2f(N);
    const UINT64 fft_nlimbs = limb_precision(fft_prec_bits);
    reps.fft_nlimbs = fft_nlimbs;
    reps.ell = ell;
    const BigComplex *powomega = fftAutoPrecomp.omega(n, fft_nlimbs);

    //compute the number of bits in the plaintext norm
    double norm_plaintext = 1;
    for (UINT64 k = 0; k < N; k++) {
        norm_plaintext += fabs(double(plaintext[k]));
    }
    reps.bits_a = UINT64(log2(norm_plaintext));

    //initialize the ciphertext
    for (UINT64 k = 0; k < 2; k++) {
        for (UINT64 i = 0; i < ell; i++) {
            for (UINT64 j = 0; j < 2; j++) {
                assert(reps.a[k][i][j] == nullptr); //trying to re-encrypt an already set TRGSW
                reps.a[k][i][j] = new_BigComplex_array(Ns2, fft_nlimbs);
            }
        }
    }

    //encrypt each component
    BigTorusParams btParams(fft_nlimbs);
    TRLweParams trLweParams(N, btParams);
    BigTorusPolynomial rescaled_plaintext(N, btParams);
    BigTorusPolynomial zero_plaintext(N, btParams);
    zero(zero_plaintext);
    TRLwe zPlusMuH_ciphertext(trLweParams);
    BigReal *zPlusMuHRj = new_BigReal_array(N, fft_nlimbs);
    for (int64_t i = 0; i < int64_t(ell); i++) { //must be signed
        for (int64_t k = 0; k < int64_t(N); k++) {
            //compute the rescaled plaintext plain / 2^32.(i+1)
            zero(rescaled_plaintext.getAT(k));
            if (i % 2 == 0) {
                rescaled_plaintext.getAT(k).limbs_end[-i / 2 - 1] = (1ul << 32ul);
            } else {
                rescaled_plaintext.getAT(k).limbs_end[-i / 2 - 1] = 1ul;
            }
            mulS64(rescaled_plaintext.getAT(k), plaintext[k], rescaled_plaintext.getAT(k), fft_nlimbs);
            //just to make sure:
            //RR::SetPrecision(fft_nlimbs*BITS_PER_LIMBS);
            //cout << "---" << endl;
            //cout << "rescPlaintext_k: " << to_RR(rescaled_plaintext.getAT(k)) << endl;
            //cout << "expected: " << to_RR(plaintext[k])/pow(to_RR(2), to_RR(32*(i+1))) << endl;
            //cout << log2Diff(to_RR(rescaled_plaintext.getAT(k)), to_RR(plaintext[k])/pow(to_RR(2), to_RR(32*(i+1)))) << endl;
        }
        for (int64_t j = 0; j < 2; j++) {
            //compute a ciphertext of zero
            native_encrypt(zPlusMuH_ciphertext, zero_plaintext, key, alpha_bits);
            //just to make sure:
            //RR::SetPrecision(fft_nlimbs*BITS_PER_LIMBS);
            //BigTorusPolynomial test(N, btParams);
            //native_phase(test, zPlusMuH_ciphertext, key, alpha_bits);
            //for (int k=0; k<int(N); k++) {
            //    cout << "test - " << log2Diff(test.getAT(i),zero_plaintext.getAT(0)) << endl;
            //}
            //translate by the rescaled plaintext
            add(zPlusMuH_ciphertext.a[j], zPlusMuH_ciphertext.a[j], rescaled_plaintext);
            /*
            //just to make sure:
            for (int64_t k = 0; k < int64_t(N); k++) {
                RR::SetPrecision(fft_nlimbs * BITS_PER_LIMBS);
                BigTorusPolynomial test(N, btParams);
                BigTorusPolynomial test2(N, btParams);
                native_phase(test, zPlusMuH_ciphertext, key, alpha_bits);
                if (j == 0)
                    fft_external_product(test2, key.key, rescaled_plaintext, 1, fft_nlimbs);
                else
                    ::copy(test2, rescaled_plaintext, fft_nlimbs);
                cout << "---" << endl;
                //cout << "phase_k: " << to_RR(test.getAT(k)) << endl;
                //cout << "expected: " << to_RR(plaintext[k]) / pow(to_RR(2), to_RR(32 * (i + 1))) << endl;
                cout << j << " " << i << " " << k << " " <<
                     log2Diff(to_RR(test.getAT(k)), (j==0?-1:1)*to_RR(test2.getAT(k))) << endl;
            }
             */

            for (int64_t l = 0; l < 2; l++) {
                //convert to BigReal
                for (int64_t k = 0; k < int64_t(N); k++) {
                    to_BigReal(zPlusMuHRj[k], zPlusMuH_ciphertext.a[l].getAT(k));
                }
                //fft to BigComplex
                iFFT(reps.a[j][i][l], zPlusMuHRj, n, powomega);
            }
        }
    }
    //free resources
    delete_BigReal_array(N, zPlusMuHRj);
}

#define DEBUG_EXTERNAL_PRODUCT 0

#if DEBUG_EXTERNAL_PRODUCT
extern int64_t *debug_plaintext;
extern TLweKey *debug_key;
#endif

void external_product(TRLwe &reps, const TRGSW &a, const TRLwe &b, int64_t out_alpha_bits) {
    const UINT64 N = a.params.N;
    assert(b.params.N == N);
    assert(reps.params.N == N);
    const UINT64 n = N * 2;
    const UINT64 Ns2 = N / 2;
    const UINT64 ell = a.ell;
    const UINT64 in_nblimbs = b.params.fixp_params.torus_limbs;
    const UINT64 out_nblimbs = limb_precision(out_alpha_bits);
    const UINT64 fft_nlimbs = a.fft_nlimbs;
    const BigComplex *powomega = fftAutoPrecomp.omega(n, fft_nlimbs);
    const BigComplex *powombar = fftAutoPrecomp.omegabar(n, fft_nlimbs);

    BigReal *poly = new_BigReal_array(N, fft_nlimbs); //decomposition in base 2^32
    BigComplex *poly_fft = new_BigComplex_array(Ns2, fft_nlimbs);

    BigTorusPolynomial tmpDec(N, b.params.fixp_params);
    BigComplex *accFFT[2];
    accFFT[0] = new_BigComplex_array(Ns2, fft_nlimbs);
    accFFT[1] = new_BigComplex_array(Ns2, fft_nlimbs);
    BigReal *acc[2];
    acc[0] = new_BigReal_array(N, fft_nlimbs);
    acc[1] = new_BigReal_array(N, fft_nlimbs);

    //set the accumulator to 0
    for (UINT64 l = 0; l < 2; l++) {
        for (UINT64 k = 0; k < Ns2; ++k) {
            zero(accFFT[l][k]);
        }
    }

#if DEBUG_EXTERNAL_PRODUCT
    //DEBUG: declare variables
    int64_t *debug_intpoly = new int64_t[N];
    TRLwe debug_acc(b.params);
    //TRLwe debug_acc2(b.params);
    BigTorusPolynomial debug_acc_phase(N, b.params.fixp_params);
    BigTorusPolynomial debug_acc2_phase(N, b.params.fixp_params);
    BigTorusPolynomial debug_hji(N, b.params.fixp_params);
    BigTorusPolynomial debug_tmp(N, b.params.fixp_params);
    zero(debug_acc);
#endif

#if DEBUG_EXTERNAL_PRODUCT
    {
        // TODO debug: at this point, verify that the phase of the accumulator
        // contains the same info as the phase of the debug accum
        native_phase(debug_acc_phase, debug_acc, *debug_key, fft_nlimbs * BITS_PER_LIMBS);
        native_phase_FFT(debug_acc2_phase, accFFT[0], accFFT[1], *debug_key, 32);
        //fft_external_product(debug_tmp, debug_plaintext, debug_acc_phase, 32, fft_nlimbs);
        //cout << "init: " << endl;
        for (int k = 0; k < int64_t(N); k++) {
            //cout << "yyyy " << k << endl;
            RR::SetPrecision(in_nblimbs * BITS_PER_LIMBS);
            RR::SetOutputPrecision(in_nblimbs * BITS_PER_LIMBS / log2(10.));
            //cout << to_RR(debug_acc_phase.getAT(k)) << endl;
            //cout << to_RR(debug_acc2_phase.getAT(k)) << endl;
            assert(log2Diff(to_RR(debug_acc2_phase.getAT(k)), to_RR(debug_acc_phase.getAT(k)))<=-1e40);
        }
    }
#endif

    for (UINT64 j = 0; j < 2; j++) {
        //tmpDec = b.a[j] + (bitDecomp_in_offset,...,bitDecomp_in_offset)
        for (UINT64 k = 0; k < N; ++k) {
            bitdecomp_signed_offset32_apply(tmpDec.getAT(k), b.a[j].getAT(k));
        }

        //decompose tmp in base 2^32
        for (int64_t i = 0; i < int64_t(ell); i++) {
            for (UINT64 k = 0; k < N; ++k) {
                int64_t polyk = bitdecomp_signed_coef32(tmpDec.getAT(k), i + 1, in_nblimbs);
#if DEBUG_EXTERNAL_PRODUCT
                debug_intpoly[k] = polyk;
#endif
                to_BigReal(poly[k], polyk, 32); //warning, in the bigReal, shift by 32
            }
#if DEBUG_EXTERNAL_PRODUCT
            // TODO debug
            zero(debug_hji);
            debug_hji.getAT(0).limbs_end[-(i + 2) / 2] = (((i % 2) == 0) ? 0x100000000ul : 0x1ul);
            fft_external_product(debug_tmp, debug_intpoly, debug_hji, 32, in_nblimbs);
            //zero(debug_acc);
            add(debug_acc.a[j], debug_acc.a[j], debug_tmp);
            //debug_acc contains debug_intpoly * hji
#endif
            //now, poly contains the decomposition coef of the plaintext over H[j][i].
            iFFT(poly_fft, poly, n, powomega); //warning, in the bigcomplex, shift by 32
            //multiply the TRLWE a.a[j][i] with poly, and add it to accum
            for (UINT64 l = 0; l < 2; l++) {
                for (UINT64 k = 0; k < Ns2; ++k) {
                    //zero(accFFT[l][k]);
                    addMulTo(accFFT[l][k], poly_fft[k], a.a[j][i][l][k]);
                }
            }
#if DEBUG_EXTERNAL_PRODUCT
            // debug: at this point, verify that the phase of the accumulator
            // contains the same info as the phase of the debug accum
            //cout << "j,i: " << j << " " << i << endl;
            native_phase_FFT(debug_acc2_phase, accFFT[0], accFFT[1], *debug_key, 32);
            native_phase(debug_acc_phase, debug_acc, *debug_key, fft_nlimbs * BITS_PER_LIMBS);
            fft_external_product(debug_tmp, debug_plaintext, debug_acc_phase, 32, fft_nlimbs);
            for (int k = 0; k < int64_t(N); k++) {
                RR::SetPrecision(in_nblimbs * BITS_PER_LIMBS);
                RR::SetOutputPrecision(in_nblimbs * BITS_PER_LIMBS / log2(10.));
                //cout << "yyyy " << k << endl;
                //cout << to_RR(debug_tmp.getAT(k)) << endl;
                //cout << to_RR(debug_acc2_phase.getAT(k)) << endl;
                assert(log2Diff(to_RR(debug_acc2_phase.getAT(k)), to_RR(debug_tmp.getAT(k))) <= -out_alpha_bits);
            }
#endif

        }
    }
#if DEBUG_EXTERNAL_PRODUCT
    //TODO debug verify that the accumulator contains b
    for (int l = 0; l < 2; l++) {
        for (int k = 0; k < int64_t(N); k++) {
            //cout << "xxxxx" << endl;
            RR::SetPrecision(in_nblimbs * BITS_PER_LIMBS);
            RR::SetOutputPrecision(in_nblimbs * BITS_PER_LIMBS / log2(10.));
            //cout << to_RR(debug_acc.a[l].getAT(k)) << endl;
            //cout << to_RR(b.a[l].getAT(k)) << endl;
            assert(log2Diff(to_RR(debug_acc.a[l].getAT(k)), to_RR(b.a[l].getAT(k))) <= -ell*32);
        }
    }
#endif

    //finally, un-FFT the accumulator into the final answer
    for (UINT64 l = 0; l < 2; l++) {
        FFT(acc[l], accFFT[l], n, powombar);
        for (UINT64 k = 0; k < N; ++k) {
            to_BigTorus(reps.a[l].getAT(k), acc[l][k], 32, out_nblimbs);
        }
    }
    //free resources
    delete_BigReal_array(N, acc[0]);
    delete_BigReal_array(N, acc[1]);
    delete_BigComplex_array(Ns2, accFFT[0]);
    delete_BigComplex_array(Ns2, accFFT[1]);
    delete_BigComplex_array(Ns2, poly_fft);
    delete_BigReal_array(N, poly);
#if DEBUG_EXTERNAL_PRODUCT
    delete[] debug_intpoly;
#endif
}


void native_phase_FFT(BigTorusPolynomial &reps, BigComplex *a, BigComplex *b, const TLweKey &key, int64_t lshift) {
    const int64_t N = reps.length;
    const int64_t n = reps.length * 2;
    const int64_t fft_nblimbs = a->real.nblimbs;
    const BigComplex *powombar = fftAutoPrecomp.omegabar(n, fft_nblimbs);
    BigReal *ra = new_BigReal_array(N, fft_nblimbs);
    BigReal *rb = new_BigReal_array(N, fft_nblimbs);
    BigTorusParams bt_params(fft_nblimbs);
    TRLweParams params(N, bt_params);
    TRLwe plain_ab(params);
    FFT(ra, a, n, powombar);
    FFT(rb, b, n, powombar);
    for (int k = 0; k < int(N); k++) {
        to_BigTorus(plain_ab.a[0].getAT(k), ra[k], lshift, fft_nblimbs);
        to_BigTorus(plain_ab.a[1].getAT(k), rb[k], lshift, fft_nblimbs);
    }
    native_phase(reps, plain_ab, key, fft_nblimbs * BITS_PER_LIMBS);
    delete_BigReal_array(N, rb);
    delete_BigReal_array(N, ra);
}

// reps = b + c.(a-b)
void cmux(TRLwe &reps, const TRGSW &c, const TRLwe &a, const TRLwe &b, int64_t out_alpha_bits) {
    TRLwe tmp(a.params);
    sub(tmp, a, b);
    external_product(tmp, c, tmp, out_alpha_bits);
    add(reps, tmp, b);
}


// reps *= X^power
void rotate_diff(TRLwe &out, const TRLwe &in, int64_t power) {
    assert(&out != &in);

    int64_t N = in.params.N;
    power = (power % (2 * N) + (2 * N)) % (2 * N);

    for (int i = 0; i < N; i++) {
        if (power + i < N) {
            copy(out.a[0].getAT(i + power), in.a[0].getAT(i));
            copy(out.a[1].getAT(i + power), in.a[1].getAT(i));
        } else if (power + i < 2 * N) {
            neg(out.a[0].getAT(i + power - N), in.a[0].getAT(i));
            neg(out.a[1].getAT(i + power - N), in.a[1].getAT(i));
        } else {
            copy(out.a[0].getAT(i + power - 2 * N), in.a[0].getAT(i));
            copy(out.a[1].getAT(i + power - 2 * N), in.a[1].getAT(i));
        }

    }
}

void rotate(TRLwe &out, const TRLwe &in, int64_t power) {

    if (&out == &in) {
        TRLwe tmp(in.params);
        copy(tmp, in);
        rotate_diff(out, tmp, power);
    } else rotate_diff(out, in, power);
}

//reps *= X^-(b-sum c_i a_i)
void blind_rotate(TRLwe &reps, int64_t b, int64_t *a, const TRGSW *c, int64_t n_in, int64_t out_alpha_bits) {

    TRLwe tmp(reps.params);
    rotate(reps, reps, -b);

    for (int i = 0; i < n_in; i++) {
        rotate(tmp, reps, a[i]);
        cmux(reps, c[i], tmp, reps, out_alpha_bits);
    }
}

TRGSW *new_TRGSW_array(UINT64 size, const TRGSWParams &params) {

    TRGSW *reps = (TRGSW *) malloc(size * sizeof(TRGSW));
    for (UINT64 i = 0; i < size; i++) {
        new(reps + i) TRGSW(params);
    }
    return reps;
}

void delete_TRGSW_array(UINT64 size, TRGSW *array) {
    for (UINT64 i = 0; i < size; i++) {
        (array + i)->~TRGSW();
    }
    free(array);
}

void fixp_internal_product(TRLwe &reps, const TRLwe &a, const TRLwe &b, const TRGSW &rk, int precision_bits) {
    const int64_t N = reps.params.N;
    const int64_t Ns2 = N / 2;
    assert_dramatically(a.params.fixp_params.level_expo > 0, "missing level info in a");
    assert_dramatically(b.params.fixp_params.level_expo > 0, "missing level info in b");
    assert_dramatically(reps.params.fixp_params.level_expo > 0, "missing level info in reps");
    assert_dramatically(a.params.N == N);
    assert_dramatically(b.params.N == N);
    assert_dramatically(rk.params.N == N);
    //compute parameters
    int64_t input_level_expo = reps.params.fixp_params.level_expo + precision_bits + 1;
    assert_dramatically(a.params.fixp_params.level_expo >= input_level_expo, "level expo of a too small");
    assert_dramatically(b.params.fixp_params.level_expo >= input_level_expo, "level expo of b too small");
    //copy and resample the input to exactly input_level
    int64_t fft_bits = 2 * input_level_expo + int64_t(log2(N)); //precision of the fft
    int64_t fft_nlimbs = limb_precision(fft_bits);
    BigReal *a1 = new_BigReal_array(N, fft_nlimbs);
    BigReal *b1 = new_BigReal_array(N, fft_nlimbs);
    BigReal *a2 = new_BigReal_array(N, fft_nlimbs);
    BigReal *b2 = new_BigReal_array(N, fft_nlimbs);
    BigReal *a1_a2 = new_BigReal_array(N, fft_nlimbs);
    BigReal *a1_b2 = new_BigReal_array(N, fft_nlimbs);
    BigReal *b1_a2 = new_BigReal_array(N, fft_nlimbs);
    BigReal *b1_b2 = new_BigReal_array(N, fft_nlimbs);
    precise_conv_toBigReal(a1, a.a[0], input_level_expo);
    precise_conv_toBigReal(b1, a.a[1], input_level_expo);
    precise_conv_toBigReal(a2, b.a[0], input_level_expo);
    precise_conv_toBigReal(b2, b.a[1], input_level_expo);
    fft_BigRealPoly_product(a1_a2, a1, a2, N, fft_nlimbs);
    fft_BigRealPoly_product(a1_b2, a1, b2, N, fft_nlimbs);
    fft_BigRealPoly_product(b1_a2, b1, a2, N, fft_nlimbs);
    fft_BigRealPoly_product(b1_b2, b1, b2, N, fft_nlimbs);

    //Now, we create the intermediate TRLWE ciphertexts
    TRLwe tmp0; //TODO find the intermediate params
    TRLwe tmp1;
    BigReal *C0 = b1_b2;
    BigReal *C1 = a1_b2;
    BigRealPoly_addTo(C1, b1_a2, N, fft_nlimbs);
    BigReal *C2 = a1_a2;
    // tmp0 = (C1 , C0) = ( a1.b2+a2.b1, b1.b2 )
    // tmp1 = (C2 , 0)  = ( a1.a2      , 0     )
    BigRealPoly_shift_toBigTorus(tmp0.a[0], C1, input_level_expo)
    BigRealPoly_shift_toBigTorus(tmp0.a[1], C0, input_level_expo)
    BigRealPoly_shift_toBigTorus(tmp1.a[0], C2, input_level_expo)
    zero(tmp1.a[1]);

    //finally return tmp0 - rk . c1
    external_product(tmp1, rk, tmp1, out_alpha_bits);
    sub(reps, tmp0, tmp1);
}
