#include <mpfr.h>
#include "BigFFT128.h"


Complex128 *precomp_iFFT128(int n, uint64_t nblimbs) {
    Complex128 *buf = (Complex128 *) malloc((n + 1) * sizeof(Complex128));
    int32_t *nn = (int32_t *) (buf);
    Complex128 *powomega = buf + 1;
    *nn = n;
    for (int i = 0; i < n; i++) {
        new(powomega + i) Complex128();
    }
    for (int i = 0; i < n; i++) {
        __float128 angle = (2. * M_PI * i) / n; //TODO
        powomega[i] = Complex128(cosf128(angle), sinf128(angle));
        //RR angle = to_RR((2.*M_PI * i) / n);
        //powomega[i] = Complex128(cos(to_RR(angle)), sin(to_RR(angle)));
    }
    return powomega;
}

void clear_precomp_iFFT128(Complex128 *powomega) {
    Complex128 *buf = powomega - 1;
    int32_t *nn = (int32_t *) buf;
    const int n = *nn;
    for (int i = 0; i < n; i++) {
        (powomega + i)->~Complex128();
    }
    free(buf);
}

Complex128 *precomp_FFT128(int n, uint64_t nblimbs) {
    Complex128 *buf = (Complex128 *) malloc((n + 1) * sizeof(Complex128));
    int32_t *nn = (int32_t *) (buf);
    Complex128 *powombar = buf + 1;
    *nn = n;
    for (int i = 0; i < n; i++) {
        new(powombar + i) Complex128();
    }
    for (int i = 0; i < n; i++) {
        __float128 angle = (M_PI * i) / n;
        powombar[i] = Complex128(cosf128(angle), -sinf128(angle));
        //RR angle = to_RR((2.*M_PI * i) / n);
        //powombar[i] = Complex128(cos(to_RR(angle)), -sin(to_RR(angle)));
    }
    return powombar;
}

void clear_precomp_FFT128(Complex128 *powombar) {
    Complex128 *buf = powombar - 1;
    int32_t *nn = (int32_t *) buf;
    const int n = *nn;
    for (int i = 0; i < n; i++) {
        (powombar + i)->~Complex128();
    }
    free(buf);
}

void iFFT(Complex128 *out, const Real128 *in, int n, const Complex128 *powomega) {
    Complex128 t1;
    Complex128 t2;
    //const int N = n/2;
    const int ns4 = n / 4;

    //interpret the input coefs as real and imaginary parts:
    //RRRRRRRRRRIIIIIIIII
    //multiply by omega^j
    for (int j = 0; j < ns4; j++)
        out[j] = Complex128(in[j], in[j + ns4]) * powomega[j];

    //at the beginning of iteration nn
    // a_{j,i} has P_{i%nn}(omega^j)
    // where j between [rev(1) and rev(3)[
    // and i between [0 and nn[
    for (int nn = ns4; nn >= 2; nn /= 2) {
        int halfnn = nn / 2;
        for (int block = 0; block < ns4; block += nn) {
            for (int off = 0; off < halfnn; off++) {
                t1 = out[block + off];
                t2 = out[block + off + halfnn];
                out[block + off] = t1 + t2;
                out[block + off + halfnn] = (t1 - t2) * powomega[(2 * (ns4 / halfnn) * off) % n];
            }
        }
    }
}

//logarithm of a power of 2
unsigned logPow2(int n) {
    return __builtin_popcount(n - 1);
}

void FFT(Real128 *out, Complex128 *in, int n, const Complex128 *powombar) {
    Complex128 t1;
    Complex128 t2;

    //const int N = n/2;
    const int ns4 = n / 4;
    //const uint64_t LOG2Ns4 = logPow2(ns4);



    //at the beginning of iteration nn
    // a_{j,i} has P_{i%nn}(omega^j)
    // where j between [rev(1) and rev(3)[
    // and i between [0 and nn[
    for (int nn = 2; nn <= ns4; nn *= 2) {
        int halfnn = nn / 2;
        for (int block = 0; block < ns4; block += nn) {
            for (int off = 0; off < halfnn; off++) {
                t1 = in[block + off];
                t2 = in[block + off + halfnn] * powombar[(2 * (ns4 / halfnn) * off) % n];
                in[block + off] = t1 + t2;
                in[block + off + halfnn] = t1 - t2;
            }
        }
    }

    //interpret the input coefs as real and imaginary parts:
    //RRRRRRRRRRIIIIIIIII
    //multiply by omega^j
    for (int j = 0; j < ns4; j++) {
        in[j] *= powombar[j];
        out[j] = in[j].real();
        out[j + ns4] = in[j].imag();
        //out[j] = in[j].real() >> LOG2Ns4;      // /ns4;  //divide by N/2
        //out[j + ns4] = in[j].imag >> LOG2Ns4;  // /ns4; //divide by N/2
    }
}
