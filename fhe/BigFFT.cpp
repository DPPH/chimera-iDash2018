#include <mpfr.h>
#include "BigFFT.h"
#include "BigReal.h"

class BigComplex {
public:
    BigReal real;
    BigReal imag;

    BigComplex(uint64_t nblimbs) : real(nblimbs), imag(nblimbs) {}
};


class BigComplexRef {
public:
    BigReal *const real;
    BigReal *const imag;

    BigComplexRef(BigReal *real, BigReal *imag) : real(real), imag(imag) {}

    BigComplexRef(const BigReal *real, const BigReal *imag) :
            real((BigReal *) real),
            imag((BigReal *) imag) {}

    BigComplexRef(BigComplex &z) :
            real(&z.real),
            imag(&z.imag) {}

    BigComplexRef(const BigComplex &z) :
            real((BigReal *) &z.real),
            imag((BigReal *) &z.imag) {}
};

void accurate_cos_sin(BigComplexRef dest, int i, int n) {
    const uint64_t nblimbs = dest.real->nblimbs;
    const uint64_t FPREC = nblimbs * BITS_PER_LIMBS;
    mpfr_t angle;
    mpfr_t rsin;
    mpfr_t rcos;
    mpfr_init2(angle, FPREC + 1);
    mpfr_init2(rsin, FPREC + 1);
    mpfr_init2(rcos, FPREC + 1);
    //compute the angle
    i = ((i % n) + n) % n;
    mpfr_const_pi(angle, MPFR_RNDN);
    mpfr_mul_si(angle, angle, 2 * i, MPFR_RNDN);
    mpfr_div_si(angle, angle, n, MPFR_RNDN);
    //compute sin and cos
    mpfr_sin_cos(rsin, rcos, angle, MPFR_RNDN);
    //shift and write the answer
    mpfr_mul_2ui(rcos, rcos, FPREC, MPFR_RNDN);
    mpfr_mul_2ui(rsin, rsin, FPREC, MPFR_RNDN);
    mpfr_get_z(dest.real->value, rcos, MPFR_RNDN);
    mpfr_get_z(dest.real->value, rsin, MPFR_RNDN);
    //cleanup
    mpfr_clear(angle);
    mpfr_clear(rsin);
    mpfr_clear(rcos);
}

void mul(BigComplexRef dest, BigComplexRef a, BigComplexRef b) {
    const uint64_t nblimbs = dest.real->nblimbs;
    BigReal t1(nblimbs);
    BigReal t2(nblimbs);
    BigReal t3(nblimbs);
    BigReal t4(nblimbs);
    mul(t1, *a.real, *b.real);
    mul(t2, *a.imag, *b.imag);
    mul(t3, *a.real, *b.imag);
    mul(t4, *a.imag, *b.real);
    sub(*dest.real, t1, t2);
    add(*dest.imag, t3, t4);
}

void mulTo(BigComplexRef dest, BigComplexRef a) {
    mul(dest, dest, a);
}

void add(BigComplexRef dest, BigComplexRef a, BigComplexRef b) {
    add(*dest.real, *a.real, *b.real);
    add(*dest.imag, *a.imag, *b.imag);
}

void sub(BigComplexRef dest, BigComplexRef a, BigComplexRef b) {
    sub(*dest.real, *a.real, *b.real);
    sub(*dest.imag, *a.imag, *b.imag);
}

void copy(BigComplexRef dest, BigComplexRef a) {
    copy(*dest.real, *a.real);
    copy(*dest.imag, *a.imag);
}

void precomp_iFFT(BigComplex *powomega, int n) {
    for (int i = 0; i < n; i++)
        accurate_cos_sin(powomega[i], i, n);
}

void precomp_FFT(BigComplex *powombar, int n) {
    for (int i = 0; i < n; i++)
        accurate_cos_sin(powombar[i], (n - i) % n, n);
}

// P -> P(omega)
// out size: n/4
// in size: n/2, //RRRRRRRRRRIIIIIIIII
void iFFT(BigComplex *out, const BigReal *in, int n, const BigComplex *powomega) {
    const uint64_t nblimbs = out[0].real.nblimbs;
    BigComplex t1(nblimbs);
    BigComplex t2(nblimbs);
    //const int N = n/2;
    const int ns4 = n / 4;

    //interpret the input coefs as real and imaginary parts:
    //RRRRRRRRRRIIIIIIIII
    //multiply by omega^j
    for (int j = 0; j < ns4; j++)
        mul(out[j], BigComplexRef(&in[j], &in[j + ns4]), powomega[j]);

    //at the beginning of iteration nn
    // a_{j,i} has P_{i%nn}(omega^j)
    // where j between [rev(1) and rev(3)[
    // and i between [0 and nn[
    for (int nn = ns4; nn >= 2; nn /= 2) {
        int halfnn = nn / 2;
        for (int block = 0; block < ns4; block += nn) {
            for (int off = 0; off < halfnn; off++) {
                copy(t1, out[block + off]);
                copy(t2, out[block + off + halfnn]);
                add(out[block + off], t1, t2);
                sub(out[block + off + halfnn], t1, t2);
                mulTo(out[block + off + halfnn], powomega[(2 * (ns4 / halfnn) * off) % n]);
            }
        }
    }
}

//logarithm of a power of 2
static unsigned logPow2(int n) {
    return __builtin_popcount(n - 1);
}

// P(omega) -> P
void FFT(BigReal *out, BigComplex *in, int n, const BigComplex *powombar) {
    const uint64_t nblimbs = out[0].nblimbs;
    BigComplex t1(nblimbs);
    BigComplex t2(nblimbs);

    //const int N = n/2;
    const int ns4 = n / 4;
    const uint64_t LOG2Ns4 = logPow2(ns4);



    //at the beginning of iteration nn
    // a_{j,i} has P_{i%nn}(omega^j)
    // where j between [rev(1) and rev(3)[
    // and i between [0 and nn[
    for (int nn = 2; nn <= ns4; nn *= 2) {
        int halfnn = nn / 2;
        for (int block = 0; block < ns4; block += nn) {
            for (int off = 0; off < halfnn; off++) {
                copy(t1, in[block + off]);
                mul(t2, in[block + off + halfnn], powombar[(2 * (ns4 / halfnn) * off) % n]);
                add(in[block + off], t1, t2);
                sub(in[block + off + halfnn], t1, t2);
            }
        }
    }

    //interpret the input coefs as real and imaginary parts:
    //RRRRRRRRRRIIIIIIIII
    //multiply by omega^j
    for (int j = 0; j < ns4; j++) {
        mulTo(in[j], powombar[j]);
        //copy(out[j], in[j].real);
        //copy(out[j + ns4], in[j].imag);
        div2ui(out[j], in[j].real, LOG2Ns4);      // /ns4;  //divide by N/2
        div2ui(out[j + ns4], in[j].imag, LOG2Ns4);  // /ns4; //divide by N/2
    }
}
