#ifndef FHE_BIGFFT_H
#define FHE_BIGFFT_H

#include "BigComplex.h"


BigComplex *precomp_iFFT(int n, uint64_t nblimbs);

void clear_precomp_iFFT(BigComplex *powomega);

BigComplex *precomp_FFT(int n, uint64_t nblimbs);

void clear_precomp_FFT(BigComplex *powombar);

// P -> P(omega)
// out size: n/4
// in size: n/2, //RRRRRRRRRRIIIIIIIII
void iFFT(BigComplex *out, const BigReal *in, int n, const BigComplex *powomega);

// P(omega) -> P
void FFT(BigReal *out, BigComplex *in, int n, const BigComplex *powombar);

#endif //FHE_BIGFFT_H
