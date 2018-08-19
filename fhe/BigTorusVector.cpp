#include "BigTorusVector.h"

using namespace std;

BigTorusVector::BigTorusVector(UINT64 length, const BigTorusParams &params) :
        btp(params),
        length(length),
        limbs_end((new UINT64[(length + 1) * params.torus_limbs]) + params.torus_limbs) {
    memset(limbs_end - params.torus_limbs, 0, (length + 1) * params.torus_limbs * sizeof(UINT64));
}

BigTorusVector::~BigTorusVector() {
    delete[] (limbs_end - btp.torus_limbs);
}

BigTorusRef BigTorusVector::getAT(UINT64 i) {
    return BigTorusRef(limbs_end + i * btp.torus_limbs, btp);
}

BigTorusRef BigTorusVector::getAT(UINT64 i) const {
    return BigTorusRef(limbs_end + i * btp.torus_limbs, btp);
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
        fixp_raw_add(
                reps.limbs_end - nreps + i * nreps,
                a.limbs_end - na + i * na,
                b.limbs_end - nb + i * nb,
                addParams);
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
        fixp_raw_sub(
                reps.limbs_end - nreps + i * nreps,
                a.limbs_end - na + i * na,
                b.limbs_end - nb + i * nb,
                addParams);
    }
    fixp_releaseAdd(addParams);
}

void subMul128(BigTorusVector &out, __int128 a, const BigTorusVector &in, const UINT64 out_limb_prec) {
    assert_dramatically(in.length == out.length, "not the good size");
    for (UINT64 i = 0; i < in.length; i++) {
        subMulS128(out.getAT(i), a, in.getAT(i), out_limb_prec);
    }
}

void subMul64(BigTorusVector &out, int64_t a, const BigTorusVector &in, const UINT64 out_limb_prec) {
    assert_dramatically(in.length == out.length, "not the good size");
    for (UINT64 i = 0; i < in.length; i++) {
        subMulS64(out.getAT(i), a, in.getAT(i), out_limb_prec);
    }
}
void copy(BigTorusVector &out, const BigTorusVector &in, UINT64 out_limbs_prec) {
    for (UINT64 i = 0; i < out.length; i++) {
        copy(out.getAT(i), in.getAT(i), out_limbs_prec);
    }
}

void random(BigTorusVector &out, UINT64 out_limbs_prec) {
    for (UINT64 i = 0; i < out.length; i++) {
        random(out.getAT(i), out_limbs_prec);
    }
}

void add_noise(BigTorusVector &out, UINT64 alpha_bits, UINT64 out_limbs_prec) {
    for (UINT64 i = 0; i < out.length; i++) {
        add_noise(out.getAT(i), alpha_bits, out_limbs_prec);
    }
}

void sub(BigTorusVector &out, const BigTorusVector &a, const BigTorusVector &b, UINT64 out_limbs_prec) {
    for (UINT64 i = 0; i < out.length; i++) {
        sub(out.getAT(i), a.getAT(i), b.getAT(i), out_limbs_prec);
    }
}

void add(BigTorusVector &out, const BigTorusVector &a, const BigTorusVector &b, UINT64 out_limbs_prec) {
    for (UINT64 i = 0; i < out.length; i++) {
        add(out.getAT(i), a.getAT(i), b.getAT(i), out_limbs_prec);
    }
}


BigTorusMatrix::BigTorusMatrix(UINT64 rows, UINT64 cols, const BigTorusParams &params) :
        rows(rows), cols(cols), params(params) {
    this->limbs_end = new UINT64[rows * cols * params.torus_limbs] + params.torus_limbs;
}

BigTorusMatrix::~BigTorusMatrix() {
    delete[] (limbs_end - params.torus_limbs);
}


void serializeBigTorusVectorContent(std::ostream &out, const BigTorusVector &value) {
    int64_t x = BIGTORUSVECTOR_SERIAL_ID;
    ostream_write_binary(out, &x, sizeof(int64_t));
    for (int i = 0; i < value.length; i++) {
        serializeBigTorusContent(out, value.getAT(i));
    }
}

void deserializeBigTorusVectorContent(std::istream &in, BigTorusVector &reps) {
    int64_t magic;

    istream_read_binary(in, &magic, sizeof(int64_t));
    assert_dramatically(magic == BIGTORUSVECTOR_SERIAL_ID);

    for (int i = 0; i < reps.length; i++) {
        deserializeBigTorusContent(in, reps.getAT(i));
    }

}


void serializeBigTorusVector(std::ostream &out, const BigTorusVector &value) {
    serializeBigTorusParams(out, value.btp);
    ostream_write_binary(out, &value.length, sizeof(UINT64));
    serializeBigTorusVectorContent(out, value);
}

std::shared_ptr<BigTorusVector> deserializeBigTorusVector(std::istream &in) {

    shared_ptr<BigTorusParams> params = deserializeBigTorusParams(in);
    store_forever(params);
    UINT64 length;
    istream_read_binary(in, &length, sizeof(UINT64));

    BigTorusVector *reps = new BigTorusVector(length, *params);
    deserializeBigTorusVectorContent(in, *reps);
    return std::shared_ptr<BigTorusVector>(reps);
}






