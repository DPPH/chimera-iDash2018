#include "TRGSW.h"
#include "BigFFT.h"

TRGSW::TRGSW(const TRGSWParams &params) : params(params) {
    for (UINT64 k = 0; k < 2; k++) {
        for (UINT64 i = 0; i < params.ell; i++) {
            for (UINT64 j = 0; j < 2; j++) {
                a[k][i][j] = new_BigComplex_array(params.N / 2, params.fixp_params.torus_limbs); //TODO
            }
        }

    }
}

TRGSW::~TRGSW() {
    for (UINT64 k = 0; k < 2; k++) {
        for (UINT64 i = 0; i < params.ell; i++) {
            for (UINT64 j = 0; j < 2; j++) {
                delete_BigComplex_array(params.N / 2, a[k][i][j]);
            }
        }
    }
}


TRGSWParams::TRGSWParams(const UINT64 N, const BigTorusParams &fixp_params) :
        TRLweParams(N, fixp_params) {}


void binary_encrypt(TRGSW &reps, const UINT64 plaintext, const TLweKey &key, UINT64 alpha_bits) {
    abort(); //TODO
}

void external_product(TRLwe &reps, TRGSW &a, TRLwe &b, UINT64 bits_a, UINT64 alpha_bits) {
    const UINT64 N = b.params.N;
    const UINT64 n = N * 2;
    const UINT64 nb_bits_decomp = alpha_bits + bits_a + 2;
    const UINT64 ell = UINT64(ceil(double(nb_bits_decomp) / double(a.params.Bgbits)));
    const UINT64 nblimbs = b.params.fixp_params.torus_limbs;

    const BigTorus bitDecomp_in_offset = ...; // sum Bg/2 Bg^i
    int64_t bitDecomp_out_offset = ...; // -Bg/2

    BigReal *poly; //decomposition in base 2^32
    const BigComplex *poly_fft;

    const BigTorusPolynomial tmp(N, b.params.fixp_params);
    const BigComplex *accFFT[2];
    const BigReal *acc[2];

    //set the accumulator to 0
    for (UINT64 l = 0; l < 2; l++) {
        for (UINT64 k = 0; k < N / 2; ++k) {
            zero(accFFT[l][k]);
        }
    }
    for (UINT64 j = 0; j < 2; j++) {
        //tmp = b.a[j] + (bitDecomp_in_offset,...,bitDecomp_in_offset)
        for (UINT64 k = 0; k < N; ++k) {
            add(tmp.getAT(k), b.a[j].getAT(k), bitDecomp_in_offset);
        }

        //decompose tmp in base 2^32
        for (UINT64 i = 0; i < ell; i++) {
            for (UINT64 k = 0; k < N; ++k) {
                int64_t polyk = bitdecomp_coef32(tmp.getAT(k), i, nblimbs) + bitDecomp_out_offset;
                to_BigReal(poly[k], polyk, 32);
            }
            iFFT(poly_fft, poly, n, powomega);
            //multiply the TRLWE a.a[j][i] with poly, and add it to accum
            for (UINT64 l = 0; l < 2; l++) {
                for (UINT64 k = 0; k < N / 2; ++k) {
                    addMulTo(accFFT[l][k], poly_fft[k], a.a[j][i][l][k]);
                }
            }
        }
    }
    //finally, un-FFT the accumulator into the final answer
    for (UINT64 l = 0; l < 2; l++) {
        FFT(acc[l], accFFT[l], n, powombar);
        for (UINT64 k = 0; k < N; ++k) {
            to_BigTorus(reps.a[l].getAT(k), acc[l][k], 32, limb_precision(alpha_bits));
        }
    }
}
