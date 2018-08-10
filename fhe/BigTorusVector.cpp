#include "BigTorusVector.h"

BigTorusVector::BigTorusVector(UINT64 length, const BigTorusParams &params) :
        btp(params),
        length(length),
        limbs(new UINT64[(length + 1) * params.torus_limbs * sizeof(UINT64)]) {

}

BigTorusVector::~BigTorusVector() {
    delete[] limbs;
}

BigTorusRef BigTorusVector::getAT(UINT64 i) {
    return BigTorusRef(limbs + i * btp.torus_limbs, &btp);
}

BigTorusRef BigTorusVector::getAT(UINT64 i) const {
    return BigTorusRef(limbs + i * btp.torus_limbs, &btp);
}

void zero(const BigTorusVector &v) {
    for (UINT64 i = 0; i < v.length; i++)
        zero(v.getAT(i));
}


void fixp_add(BigTorusVector &reps, const BigTorusVector &a, const BigTorusVector &b, UINT64 out_precision_bits) {
    const int64_t nreps = reps.btp.torus_limbs;
    const int64_t na = a.btp.torus_limbs;
    const int64_t nb = b.btp.torus_limbs;
    assert_dramatically(reps.length == a.length && reps.length == b.length, "wrong dimensions");
    fixp_add_params addParams;
    if (out_precision_bits == NA) out_precision_bits = reps.btp.torus_limbs * BITS_PER_LIMBS;
    fixp_prepareAdd(addParams, reps.btp, a.btp, b.btp, out_precision_bits);
    for (UINT64 i = 0; i < a.length; i++) {
        fixp_raw_add(reps.limbs + i * nreps, a.limbs + i * na, b.limbs + i * nb, addParams);
    }
    fixp_releaseAdd(addParams);
}

void fixp_sub(BigTorusVector &reps, const BigTorusVector &a, const BigTorusVector &b, UINT64 out_precision_bits) {
    const int64_t nreps = reps.btp.torus_limbs;
    const int64_t na = a.btp.torus_limbs;
    const int64_t nb = b.btp.torus_limbs;
    assert_dramatically(reps.length == a.length && reps.length == b.length, "wrong dimensions");
    fixp_add_params addParams;
    if (out_precision_bits == NA) out_precision_bits = reps.btp.torus_limbs * BITS_PER_LIMBS;
    fixp_prepareAdd(addParams, reps.btp, a.btp, b.btp, out_precision_bits);
    for (UINT64 i = 0; i < a.length; i++) {
        fixp_raw_sub(reps.limbs + i * nreps, a.limbs + i * na, b.limbs + i * nb, addParams);
    }
    fixp_releaseAdd(addParams);
}

void subMul(BigTorusVector &out, __int128 a, const BigTorusVector &in, const UINT64 out_limb_prec) {
    assert_dramatically(in.length == out.length, "not the good size");
    for (int i = 0; i < in.length; i++) {
        subMul(out.getAT(i), a, in.getAT(i), out_limb_prec);
    }

}


BigTorusMatrix::BigTorusMatrix(UINT64 rows, UINT64 cols, const BigTorusParams *params) :
        rows(rows), cols(cols), params(params) {
    this->limbs = new UINT64[rows * cols * params->torus_limbs];
}

BigTorusMatrix::~BigTorusMatrix() {
    delete[] limbs;
}
