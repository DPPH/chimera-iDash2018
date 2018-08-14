#include "TRGSW.h"
#include "BigFFT.h"

TRGSW::TRGSW(const TRGSWParams &params) : params(params) {

    for (UINT64 i = 0; i < params.ell; i++) {
        for (UINT64 j = 0; i < 2; i++) {

            a[i][j] = new_BigComplex_array(params.N / 2, params.fixp_params.torus_limbs); //TODO
        }
    }

}

TRGSW::~TRGSW() {
    for (UINT64 i = 0; i < params.ell; i++) {
        for (UINT64 j = 0; i < 2; i++) {

            delete_BigComplex_array(params.N / 2, a[i][j]);
        }
    }
}


TRGSWParams::TRGSWParams(const UINT64 N, const BigTorusParams &fixp_params) : TRLweParams(N, fixp_params) {

}


void binary_encrypt(TRGSW &reps, const UINT64 plaintext, const TLweKey &key, UINT64 alpha_bits) {
    abort();
}

void external_product(TRLwe &reps, TRGSW &a, TRLwe &b, UINT64 bits_a, UINT64 alpha_bits) {
    UINT64 n = b.params.N * 2;
    UINT64 nb_bits_decomp = alpha_bits + bits_a;
    UINT64 ell = ceil(double(nb_bits_decomp) / double(a.params.Bgbits));
    BigReal *poly;
    const BigTorus &bitDecomp_in_offset = bitDecomp_in_offset; // sum Bg/2 Bg^i
    int64_t bitDecomp_out_offset = bitDecomp_out_offset; // -Bg/2
    const BigTorusPolynomial tmp(b.params.N, b.params.fixp_params);
    const BigComplex *acc[2];
    const BigComplex *poly_fft;

    for (UINT64 j = 0; j < 2; j++) {
        for (UINT64 k = 0; k < b.params.N; ++k) {
            add(tmp.getAT(k), b.a[j].getAT(k), bitDecomp_in_offset);
        }

        for (UINT64 i = 0; i < ell; i++) {
            for (UINT64 k = 0; k < b.params.N; ++k) {
                int polyk = bitdecomp_coef32(tmp.getAT(k), i, b.params.fixp_params.torus_limbs);
                to_BigReal(poly[k], polyk, 32);
            }
            iFFT(poly_fft, poly, n, powomega);
            addmulto(acc, poly_fft, a.a[i]);
        }

    }




}
