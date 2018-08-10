#ifndef FHE_BIGFFT128_H
#define FHE_BIGFFT128_H

#include <complex>
#include "commons.h"

typedef __float128 Real128;
typedef std::complex<Real128> Complex128;


/**
 * precompute an iFFT structure for logical n (=2*N if we are mod X^N+1) and limb precision nblimbs
 * Although this structure is mostly opaque, in this poc, the returned array contains all complex powers of unity
 * TODO!
 */
Complex128 *precomp_iFFT128(int n, UINT64 nblimbs);

/** releases the precomputed an iFFT structure */
void clear_precomp_iFFT128(Complex128 *powomega);

/**
 * precompute an FFT structure for logical n (=2*N if we are mod X^N+1) and limb precision nblimbs
 * Although this structure is mostly opaque, in this poc, the returned array contains all complex powers of unity
 * in reverse order
 * TODO!
 */
Complex128 *precomp_FFT128(int n, UINT64 nblimbs);

/** releases the precomputed an iFFT structure */
void clear_precomp_FFT128(Complex128 *powombar);

/**
 * computes the iFFT of a real polynomial P(X) mod X^N+1
 * that is, the complex array [P(omega_1),...,P(omega_N/2)]
 * WARNING: in order to ensure that no overflow occurs, all input coefficients should be < n/4 in absolute value
 * @param out the complex array (of size N/2 = n/4)
 * @param in the real coefficients of P (of size N = n/2)
 * @param n the logical dimension (= 2N)
 * @param powomega the precomputed iFFT structure
 */
void iFFT(Complex128 *out, const Real128 *in, int n, const Complex128 *powomega);

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
void FFT(Real128 *out, Complex128 *in, int n, const Complex128 *powombar);

#endif //FHE_BIGFFT128_H
