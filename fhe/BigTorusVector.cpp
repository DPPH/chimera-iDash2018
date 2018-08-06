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


