#include <cassert>
#include "TRGSW.h"
#include "BigFFT.h"

TRGSW::TRGSW(const TRGSWParams &params) : params(params), bits_a(0), fft_nlimbs(0), ell(0) {
    for (UINT64 k = 0; k < 2; k++) {
        for (UINT64 i = 0; i < params.max_ell; i++) {
            for (UINT64 j = 0; j < 2; j++) {
                a[k][i][j] = new_BigComplex_array(params.N / 2, params.fixp_params.torus_limbs); //TODO
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
    const UINT64 nb_bits_decomp = alpha_bits + 2;
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
    BigTorusPolynomial rescaled_plaintext(N, fft_nlimbs);
    BigTorusPolynomial zero_plaintext(N, fft_nlimbs);
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
        }
        for (int64_t j = 0; j < 2; j++) {
            //compute a ciphertext of zero
            native_encrypt(zPlusMuH_ciphertext, zero_plaintext, key, alpha_bits);
            //translate by the rescaled plaintext
            add(zPlusMuH_ciphertext.a[j], zPlusMuH_ciphertext.a[j], rescaled_plaintext);

            for (int64_t l = 0; l < 2; l++) {
                //convert to BigReal
                for (int64_t k = 0; k < int64_t(N); k++) {
                    to_BigReal(zPlusMuHRj[k], zPlusMuH_ciphertext.a[j].getAT(k));
                }
                //fft to BigComplex
                iFFT(reps.a[j][i][l], zPlusMuHRj, n, powomega);
            }
        }
    }
    //free resources
    delete_BigReal_array(N, zPlusMuHRj);
}

void external_product(TRLwe &reps, TRGSW &a, TRLwe &b, UINT64 out_alpha_bits) {
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
    BigComplex *accFFT[2] = {new_BigComplex_array(Ns2, fft_nlimbs), new_BigComplex_array(Ns2, fft_nlimbs)};
    BigReal *acc[2] = {new_BigReal_array(N, fft_nlimbs), new_BigReal_array(N, fft_nlimbs)};;

    //set the accumulator to 0
    for (UINT64 l = 0; l < 2; l++) {
        for (UINT64 k = 0; k < Ns2; ++k) {
            zero(accFFT[l][k]);
        }
    }
    for (UINT64 j = 0; j < 2; j++) {
        //tmp = b.a[j] + (bitDecomp_in_offset,...,bitDecomp_in_offset)
        for (UINT64 k = 0; k < N; ++k) {
            bitdecomp_signed_offset32_apply(tmpDec.getAT(k), b.a[j].getAT(k));
        }

        //decompose tmp in base 2^32
        for (UINT64 i = 0; i < ell; i++) {
            for (UINT64 k = 0; k < N; ++k) {
                int64_t polyk = bitdecomp_signed_coef32(tmpDec.getAT(k), i, in_nblimbs);
                to_BigReal(poly[k], polyk, 32);
            }
            iFFT(poly_fft, poly, n, powomega);
            //multiply the TRLWE a.a[j][i] with poly, and add it to accum
            for (UINT64 l = 0; l < 2; l++) {
                for (UINT64 k = 0; k < Ns2; ++k) {
                    addMulTo(accFFT[l][k], poly_fft[k], a.a[j][i][l][k]);
                }
            }
        }
    }
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
}
