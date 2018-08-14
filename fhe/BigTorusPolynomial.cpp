#include "BigTorusPolynomial.h"
#include "BigComplex.h"

BigTorusPolynomial::BigTorusPolynomial(UINT64 N, const BigTorusParams &params) : BigTorusVector(N, params) {
}

BigTorusPolynomial::~BigTorusPolynomial() {

}

void const_poly(BigTorusPolynomial &out, const BigTorusRef &in, UINT64 out_limb_prec) {
    copy(out.getAT(0), in, out_limb_prec);
    for (UINT64 i = 1; i < out.length; i++) {
        zero(out.getAT(i));
    }
}

void naive_external_product(BigTorusPolynomial &out, int64_t *a, const BigTorusPolynomial &b, UINT64 out_limb_prec) {
    const UINT64 N = out.length;
    if (out_limb_prec == NA) out_limb_prec = out.btp.torus_limbs;
    //in this function, input precision is min(out_limb_prec + 1, in_limb_prec)
    assert_dramatically(out_limb_prec <= out.btp.torus_limbs, "output precision requested too high");
    assert_dramatically(out_limb_prec <= b.btp.torus_limbs, "input precision too small vrt. output prec");
    const UINT64 in_limb_prec = std::min(b.btp.torus_limbs, out_limb_prec + 1);

    BigTorus ri(b.btp);
    BigTorus tmp(b.btp);
    for (UINT64 i = 0; i < N; i++) {
        zero(ri);
        for (UINT64 j = 0; j <= i; j++) {
            //ri += poly1[j]*poly2[i-j];
            mulS64(tmp, a[j], b.getAT(i - j), in_limb_prec);
            add(ri, ri, tmp, in_limb_prec);
        }
        for (UINT64 j = i + 1; j < N; j++) {
            //ri -= poly1[j]*poly2[N+i-j];
            mulS64(tmp, a[j], b.getAT(N + i - j), in_limb_prec);
            sub(ri, ri, tmp, in_limb_prec);
        }
        copy(out.getAT(i), ri, out_limb_prec);
    }
}

void fft_internal_product(BigTorusPolynomial &out, const BigTorusPolynomial &a, const BigTorusPolynomial &b,
                          UINT64 out_limb_prec) {

    assert_dramatically(a.length == b.length, "not the same size of a and b");
    UINT64 N = b.length;
    UINT64 n = 2 * b.length;
    UINT64 prec_bits = b.btp.torus_limbs * BITS_PER_LIMBS + a.btp.torus_limbs * BITS_PER_LIMBS + (UINT64) log2f(N);
    UINT64 nblimbs = limb_precision(prec_bits);


    BigReal *rb = new_BigReal_array(N, nblimbs);
    BigReal *ra = new_BigReal_array(N, nblimbs);
    BigComplex *cb = new_BigComplex_array(N / 2, nblimbs);
    BigComplex *ca = new_BigComplex_array(N / 2, nblimbs);

    for (UINT64 i = 0; i < N; i++) {
        to_BigReal(rb[i], b.getAT(i));
        to_BigReal(ra[i], a.getAT(i));

    }

    iFFT(cb, rb, n, powomega);
    iFFT(ca, ra, n, powomega);
    for (UINT64 i = 0; i < N / 2; i++) {
        mulTo(cb[i], ca[i]);
    }

    FFT(rb, cb, n, powombar);
    for (UINT64 i = 0; i < N; i++) {
        to_BigTorus(out.getAT(i), rb[i], 0, out_limb_prec); //TODO
    }

}

void fft_semi_external_product(BigTorusPolynomial &out, const BigComplex *ca, const BigTorusPolynomial &b,
                               const UINT64 bits_a, UINT64 out_limb_prec) {
    UINT64 N = b.length;
    UINT64 n = 2 * b.length;
    UINT64 nblimbs = ca->real.nblimbs;


    BigReal *rb = new_BigReal_array(N, nblimbs);
    BigComplex *cb = new_BigComplex_array(N / 2, nblimbs);

    for (UINT64 i = 0; i < N; i++) {
        to_BigReal(rb[i], b.getAT(i));

    }

    iFFT(cb, rb, n, powomega);

    for (UINT64 i = 0; i < N / 2; i++) {
        mulTo(cb[i], ca[i]);
    }

    FFT(rb, cb, n, powombar);
    for (UINT64 i = 0; i < N; i++) {
        to_BigTorus(out.getAT(i), rb[i], bits_a, out_limb_prec); //TODO
    }
}

void fft_external_product(BigTorusPolynomial &out, int64_t *a, const BigTorusPolynomial &b, const UINT64 bits_a,
                          UINT64 out_limb_prec) {

    UINT64 N = b.length;
    UINT64 n = 2 * b.length;
    UINT64 prec_bits = b.btp.torus_limbs * BITS_PER_LIMBS + bits_a + (UINT64) log2f(N);
    UINT64 nblimbs = limb_precision(prec_bits);


    BigReal *rb = new_BigReal_array(N, nblimbs);
    BigReal *ra = new_BigReal_array(N, nblimbs);
    BigComplex *cb = new_BigComplex_array(N / 2, nblimbs);
    BigComplex *ca = new_BigComplex_array(N / 2, nblimbs);

    for (UINT64 i = 0; i < N; i++) {
        to_BigReal(rb[i], b.getAT(i));
        to_BigReal(ra[i], a[i], bits_a);

    }

    iFFT(cb, rb, n, powomega);
    iFFT(ca, ra, n, powomega);
    for (UINT64 i = 0; i < N / 2; i++) {
        mulTo(cb[i], ca[i]);
    }

    FFT(rb, cb, n, powombar);
    for (UINT64 i = 0; i < N; i++) {
        to_BigTorus(out.getAT(i), rb[i], bits_a, out_limb_prec); //TODO
    }

}
