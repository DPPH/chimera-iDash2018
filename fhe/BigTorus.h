#ifndef FHE_BIGTORUS_H
#define FHE_BIGTORUS_H

#include <cstdint>

class BigTorusParams {
public:
    uint64_t torus_limbs; ///< total number of limbs

    BigTorusParams(uint64_t torus_limbs);
};

static_assert(sizeof(BigTorusParams) == 8, "Bad compiler");

class BigTorus {
public:
    uint64_t *limbs_raw;
    const BigTorusParams *const params;

    BigTorus(const BigTorusParams *params) :
            params(params) {
        limbs_raw = new uint64_t[params->torus_limbs];
    }

    ~BigTorus() {
        delete[] limbs_raw;
    }
};

class BigTorusVector {
public:
    uint64_t length;
    uint64_t *limbs_raw;
    BigTorusParams *params;
};

class BigTorusMatrix {
public:
    uint64_t rows;
    uint64_t cols;
    uint64_t *limbs_raw;
    BigTorusParams *params;
};

class BigTorusRef {
public:
    uint64_t *limbs_raw;
    const BigTorusParams *const params;

    BigTorusRef(uint64_t *limbs_raw, const BigTorusParams *params);

    BigTorusRef(const BigTorus &torus);

    BigTorusRef(BigTorus &torus);
};

void bigTorusRawScale(uint64_t *limbs, int64_t coef, uint64_t nblimbs);

void bigTorusScale(const BigTorusRef &x, int64_t coef);


#endif //FHE_BIGTORUS_H
