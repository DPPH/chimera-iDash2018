#ifndef FHE_BIGTORUSPOLYNOMIAL_H
#define FHE_BIGTORUSPOLYNOMIAL_H


#include "BigTorusVector.h"

class BigTorusPolynomial : public BigTorusVector {
public:
    BigTorusPolynomial(UINT64 N, const BigTorusParams &params);

    ~BigTorusPolynomial();
};

/** @brief generate a constant polynomial equal to in */
void const_poly(BigTorusPolynomial &out, const BigTorusRef &in, UINT64 out_limb_prec);

void naive_external_product(BigTorusPolynomial &out, int64_t *a, const BigTorusPolynomial &b, UINT64 out_limb_prec) {
    const UINT64 N = out.length;
    if (out_limb_prec == NA) out_limb_prec = out.btp.torus_limbs;
    //in this function, input precision is min(out_limb_prec + 1, in_limb_prec)
    assert_dramatically(out_limb_prec <= out.btp.torus_limbs, "output precision requested too high");
    assert_dramatically(out_limb_prec <= b.btp.torus_limbs, "input precision too small vrt. output prec");
    const UINT64 in_limb_prec = std::min(b.btp.torus_limbs, out_limb_prec + 1);

    BigTorus ri(&b.btp);
    BigTorus tmp(&b.btp);
    for (UINT64 i = 0; i < N; i++) {
        zero(ri);
        for (UINT64 j = 0; j <= i; j++) {
            //ri += poly1[j]*poly2[i-j];
            mul64(tmp, a[j], b.getAT(i - j), in_limb_prec);
            add(ri, ri, tmp, in_limb_prec);
        }
        for (UINT64 j = i + 1; j < N; j++) {
            //ri -= poly1[j]*poly2[N+i-j];
            mul64(tmp, a[j], b.getAT(N + i - j), in_limb_prec);
            sub(ri, ri, tmp, in_limb_prec);
        }
        copy(out.getAT(i), ri, out_limb_prec);
    }
}

#endif //FHE_BIGTORUSPOLYNOMIAL_H
