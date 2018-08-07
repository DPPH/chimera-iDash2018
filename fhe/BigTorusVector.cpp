#include "BigTorusVector.h"

BigTorusVector::BigTorusVector(uint64_t length, const BigTorusParams &params) :
        btp(params),
        length(length),
        limbs(new uint64_t[(length + 1) * params.torus_limbs * sizeof(uint64_t)]) {

}

BigTorusVector::~BigTorusVector() {
    delete[] limbs;
}

BigTorusRef BigTorusVector::getAT(uint64_t i) {
    return BigTorusRef(limbs + i * btp.torus_limbs, &btp);
}

BigTorusRef BigTorusVector::getAT(uint64_t i) const {
    return BigTorusRef(limbs + i * btp.torus_limbs, &btp);
}

void zero(const BigTorusVector &v) {
    for (uint64_t i = 0; i < v.length; i++)
        zero(v.getAT(i));
}


void fixp_add(BigTorusVector &reps, const BigTorusVector &a, const BigTorusVector &b, uint64_t out_precision_bits) {
    const int64_t nreps = reps.btp.torus_limbs;
    const int64_t na = a.btp.torus_limbs;
    const int64_t nb = b.btp.torus_limbs;
    assert_dramatically(reps.length == a.length && reps.length == b.length, "wrong dimensions");
    fixp_add_params addParams;
    if (out_precision_bits == NA) out_precision_bits = reps.btp.torus_limbs * BITS_PER_LIMBS;
    fixp_prepareAdd(addParams, reps.btp, a.btp, b.btp, out_precision_bits);
    for (uint64_t i = 0; i < a.length; i++) {
        fixp_raw_add(reps.limbs + i * nreps, a.limbs + i * na, b.limbs + i * nb, addParams);
    }
    fixp_releaseAdd(addParams);
}

void fixp_sub(BigTorusVector &reps, const BigTorusVector &a, const BigTorusVector &b, uint64_t out_precision_bits) {
    const int64_t nreps = reps.btp.torus_limbs;
    const int64_t na = a.btp.torus_limbs;
    const int64_t nb = b.btp.torus_limbs;
    assert_dramatically(reps.length == a.length && reps.length == b.length, "wrong dimensions");
    fixp_add_params addParams;
    if (out_precision_bits == NA) out_precision_bits = reps.btp.torus_limbs * BITS_PER_LIMBS;
    fixp_prepareAdd(addParams, reps.btp, a.btp, b.btp, out_precision_bits);
    for (uint64_t i = 0; i < a.length; i++) {
        fixp_raw_sub(reps.limbs + i * nreps, a.limbs + i * na, b.limbs + i * nb, addParams);
    }
    fixp_releaseAdd(addParams);
}


BigTorusMatrix::BigTorusMatrix(uint64_t rows, uint64_t cols, const BigTorusParams *params) :
        rows(rows), cols(cols), params(params) {
    this->limbs = new uint64_t[rows * cols * params->torus_limbs];
}

BigTorusMatrix::~BigTorusMatrix() {
    delete[] limbs;
}
