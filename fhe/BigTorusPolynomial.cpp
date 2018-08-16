#include <assert.h>
#include "BigTorusPolynomial.h"
#include "BigComplex.h"
#include "BigFFT.h"

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
    BigComplex *powomega = fftAutoPrecomp.omega(n, nblimbs);
    BigComplex *powombar = fftAutoPrecomp.omegabar(n, nblimbs);
    iFFT(cb, rb, n, powomega);
    iFFT(ca, ra, n, powomega);
    for (UINT64 i = 0; i < N / 2; i++) {
        mulTo(cb[i], ca[i]);
    }

    FFT(rb, cb, n, powombar);
    for (UINT64 i = 0; i < N; i++) {
        to_BigTorus(out.getAT(i), rb[i], 0, out_limb_prec); //TODO
    }

    delete_BigReal_array(N, rb);
    delete_BigReal_array(N, ra);
    delete_BigComplex_array(N / 2, cb);
    delete_BigComplex_array(N / 2, ca);
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
    BigComplex *powomega = fftAutoPrecomp.omega(n, nblimbs);
    BigComplex *powombar = fftAutoPrecomp.omegabar(n, nblimbs);
    iFFT(cb, rb, n, powomega);

    for (UINT64 i = 0; i < N / 2; i++) {
        mulTo(cb[i], ca[i]);
    }

    FFT(rb, cb, n, powombar);
    for (UINT64 i = 0; i < N; i++) {
        to_BigTorus(out.getAT(i), rb[i], bits_a, out_limb_prec); //TODO
    }
}

using namespace NTL; //todo remove
using namespace std;

double log2Diff(const RR &a, const RR &b);

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
    BigComplex *powomega = fftAutoPrecomp.omega(n, nblimbs);
    BigComplex *powombar = fftAutoPrecomp.omegabar(n, nblimbs);

    iFFT(cb, rb, n, powomega);
    iFFT(ca, ra, n, powomega);
    for (UINT64 i = 0; i < N / 2; i++) {
        mulTo(cb[i], ca[i]);
    }

    FFT(rb, cb, n, powombar);
    for (UINT64 i = 0; i < N; i++) {
        to_BigTorus(out.getAT(i), rb[i], bits_a, out_limb_prec);
    }

    delete_BigComplex_array(N / 2, cb);
    delete_BigComplex_array(N / 2, ca);
    delete_BigReal_array(N, rb);
    delete_BigReal_array(N, ra);
}

void add(BigTorusPolynomial &reps, const BigTorusPolynomial &a, const BigTorusPolynomial &b) {
    const UINT64 N = reps.length;
    assert(a.length == N);
    assert(b.length == N);
    for (UINT64 i = 0; i < N; i++) {
        add(reps.getAT(i), a.getAT(i), b.getAT(i));
    }
}

void sub(BigTorusPolynomial &reps, const BigTorusPolynomial &a, const BigTorusPolynomial &b) {
    const UINT64 N = reps.length;
    assert(a.length == N);
    assert(b.length == N);
    for (UINT64 i = 0; i < N; i++) {
        sub(reps.getAT(i), a.getAT(i), b.getAT(i));
    }
}

void neg(BigTorusPolynomial &out, const BigTorusPolynomial &in) {

    assert(in.length == out.length);
    for (UINT64 i = 0; i < in.length; i++) {
        neg(out.getAT(i), in.getAT(i));
    }
}

void copy(BigTorusPolynomial &out, const BigTorusPolynomial &in) {
    assert(in.length == out.length);
    for (UINT64 i = 0; i < in.length; i++) {
        copy(out.getAT(i), in.getAT(i));
    }
}
