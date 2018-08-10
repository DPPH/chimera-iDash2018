#include <mpfr.h>
#include "BigFFT.h"
#include "BigReal.h"
#include "BigComplex.h"


BigComplex *precomp_iFFT(int n, UINT64 nblimbs) {
    BigComplex *buf = (BigComplex *) malloc((n + 1) * sizeof(BigComplex));
    int32_t *nn = (int32_t *) (buf);
    BigComplex *powomega = buf + 1;
    *nn = n;
    for (int i = 0; i < n; i++) {
        new(powomega + i) BigComplex(nblimbs);
    }
    for (int i = 0; i < n; i++) {
        accurate_power_unity(powomega[i], i, n);
    }
    return powomega;
}

void clear_precomp_iFFT(BigComplex *powomega) {
    BigComplex *buf = powomega - 1;
    int32_t *nn = (int32_t *) buf;
    const int n = *nn;
    for (int i = 0; i < n; i++) {
        (powomega + i)->~BigComplex();
    }
    free(buf);
}

BigComplex *precomp_FFT(int n, UINT64 nblimbs) {
    BigComplex *buf = (BigComplex *) malloc((n + 1) * sizeof(BigComplex));
    int32_t *nn = (int32_t *) (buf);
    BigComplex *powombar = buf + 1;
    *nn = n;
    for (int i = 0; i < n; i++) {
        new(powombar + i) BigComplex(nblimbs);
    }
    for (int i = 0; i < n; i++)
        accurate_power_unity(powombar[i], (n - i) % n, n);
    return powombar;
}

void clear_precomp_FFT(BigComplex *powombar) {
    BigComplex *buf = powombar - 1;
    int32_t *nn = (int32_t *) buf;
    const int n = *nn;
    for (int i = 0; i < n; i++) {
        (powombar + i)->~BigComplex();
    }
    free(buf);
}

void iFFT(BigComplex *out, const BigReal *in, int n, const BigComplex *powomega) {
    const UINT64 nblimbs = out[0].real.nblimbs;
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

void FFT(BigReal *out, BigComplex *in, int n, const BigComplex *powombar) {
    const UINT64 nblimbs = out[0].nblimbs;
    BigComplex t1(nblimbs);
    BigComplex t2(nblimbs);

    //const int N = n/2;
    const int ns4 = n / 4;
    const UINT64 LOG2Ns4 = logPow2(ns4);



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
