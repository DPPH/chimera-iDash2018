#include "BigTorusPolynomial.h"

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
