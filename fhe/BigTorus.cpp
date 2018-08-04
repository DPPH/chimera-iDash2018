#include <gmp.h>
#include "BigTorus.h"
#include "commons.h"

BigTorusParams::BigTorusParams(uint64_t torus_limbs) : torus_limbs(torus_limbs) {}

BigTorusRef::BigTorusRef(uint64_t *limbs_raw, const BigTorusParams *params) : limbs_raw(limbs_raw), params(params) {}

BigTorusRef::BigTorusRef(const BigTorus &torus) : limbs_raw(torus.limbs_raw), params(torus.params) {}

BigTorusRef::BigTorusRef(BigTorus &torus) : limbs_raw(torus.limbs_raw), params(torus.params) {}

void bigTorusRawScale(uint64_t *limbs, int64_t coef, uint64_t nblimbs) {
    assert_dramatically(coef > 0, "negative scale not supported");
    mpn_mul_1(limbs, limbs, nblimbs, coef);
}

void bigTorusScale(const BigTorusRef &x, int64_t coef) {
    bigTorusRawScale(x.limbs_raw, coef, x.params->torus_limbs);
}
