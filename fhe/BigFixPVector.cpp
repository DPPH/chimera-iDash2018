#include "BigFixPVector.h"

BigFixPVector::BigFixPVector(uint64_t length, const BigFixPParams &params) :
        BigTorusVector(length, params), bfp(params) {
    assert(&btp == &bfp); //btp anf bfp are the same object
}

BigFixPVector::~BigFixPVector() {
}

BigFixPRef BigFixPVector::getAF(uint64_t i) {
    return BigFixPRef(limbs + i * bfp.torus_limbs, &bfp);
}

BigFixPRef BigFixPVector::getAF(uint64_t i) const {
    return BigFixPRef(limbs + i * bfp.torus_limbs, &bfp);
}

BigFixPMatrix::BigFixPMatrix(uint64_t rows, uint64_t cols, const BigFixPParams *params) : rows(rows), cols(cols),
                                                                                          params(params) {
    limbs_raw = new uint64_t[rows * cols * params->torus_limbs];
}

BigFixPMatrix::~BigFixPMatrix() {
    delete[] limbs_raw;
}

BigFixPRef BigFixPMatrix::operator()(uint64_t i, uint64_t j) {
    return BigFixPRef(limbs_raw + (i * cols + j) * params->torus_limbs, params);
}

void clear(BigFixPVector &reps) {
    const int64_t nblimbs = reps.bfp.torus_limbs;
    for (uint64_t i = 0; i < reps.length; i++) {
        fixPRawClear(reps.limbs + i * nblimbs, nblimbs);
    }
}

void add(BigFixPVector &reps, const BigFixPVector &a, const BigFixPVector &b, uint64_t out_precision_bits) {
    const int64_t nreps = reps.bfp.torus_limbs;
    const int64_t na = a.bfp.torus_limbs;
    const int64_t nb = b.bfp.torus_limbs;
    assert_dramatically(reps.length == a.length && reps.length == b.length, "wrong dimensions");
    BigFixPAddParams addParams;
    if (out_precision_bits == NA) out_precision_bits = reps.bfp.torus_limbs * BITS_PER_LIMBS;
    prepareAdd(addParams, reps.bfp, a.bfp, b.bfp, out_precision_bits);
    for (uint64_t i = 0; i < a.length; i++) {
        fixPRawAdd(reps.limbs + i * nreps, a.limbs + i * na, b.limbs + i * nb, addParams);
    }
    releaseAdd(addParams);
}

void sub(BigFixPVector &reps, const BigFixPVector &a, const BigFixPVector &b, uint64_t out_precision_bits) {
    const int64_t nreps = reps.bfp.torus_limbs;
    const int64_t na = a.bfp.torus_limbs;
    const int64_t nb = b.bfp.torus_limbs;
    assert_dramatically(reps.length == a.length && reps.length == b.length, "wrong dimensions");
    BigFixPAddParams addParams;
    if (out_precision_bits == NA) out_precision_bits = reps.bfp.torus_limbs * BITS_PER_LIMBS;
    prepareAdd(addParams, reps.bfp, a.bfp, b.bfp, out_precision_bits);
    for (uint64_t i = 0; i < a.length; i++) {
        fixPRawSub(reps.limbs + i * nreps, a.limbs + i * na, b.limbs + i * nb, addParams);
    }
    releaseAdd(addParams);
}
