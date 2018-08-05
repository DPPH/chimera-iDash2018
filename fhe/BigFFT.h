#ifndef FHE_BIGFFT_H
#define FHE_BIGFFT_H

#include "BigComplex.h"


/**
 * precompute an iFFT structure for logical n (=2*N if we are mod X^N+1) and limb precision nblimbs
 * Although this structure is mostly opaque, in this poc, the returned array contains all complex powers of unity
 */
BigComplex *precomp_iFFT(int n, uint64_t nblimbs);

/** releases the precomputed an iFFT structure */
void clear_precomp_iFFT(BigComplex *powomega);

/**
 * precompute an FFT structure for logical n (=2*N if we are mod X^N+1) and limb precision nblimbs
 * Although this structure is mostly opaque, in this poc, the returned array contains all complex powers of unity
 * in reverse order
 */
BigComplex *precomp_FFT(int n, uint64_t nblimbs);

/** releases the precomputed an iFFT structure */
void clear_precomp_FFT(BigComplex *powombar);

/**
 * computes the iFFT of a real polynomial P(X) mod X^N+1
 * that is, the complex array [P(omega_1),...,P(omega_N/2)]
 * WARNING: in order to ensure that no overflow occurs, all input coefficients should be < n/4 in absolute value
 * @param out the complex array (of size N/2 = n/4)
 * @param in the real coefficients of P (of size N = n/2)
 * @param n the logical dimension (= 2N)
 * @param powomega the precomputed iFFT structure
 */
void iFFT(BigComplex *out, const BigReal *in, int n, const BigComplex *powomega);

/**
 * computes the FFT of a complex array of N/2 values,
 * that is, the unique polynomial P in R[X]/X^N
 * that is, the complex array [P(omega_1),...,P(omega_N/2)]
 * WARNING: in order to ensure that no overflow occurs, all input coefficients should be < n/4 in absolute value
 * @param out the real coefficients of P (of size N = n/2)
 * @param in the complex array (of size N/2 = n/4)
 * @param n the logical dimension (= 2N)
 * @param powombar the precomputed FFT structure
 */
void FFT(BigReal *out, BigComplex *in, int n, const BigComplex *powombar);

#endif //FHE_BIGFFT_H
