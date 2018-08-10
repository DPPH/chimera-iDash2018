#ifndef FHE_ARITHMETIC_H
#define FHE_ARITHMETIC_H

#include "BigTorusVector.h"
#include "BigTorus.h"

void to_fixP(BigTorusRef reps, const NTL::RR &a);

NTL::RR fixp_to_RR(const BigTorusRef &a);

NTL::RR to_RR(const BigTorusRef &a);

void fixp_tAb_prod(BigTorusVector &res, BigTorusMatrix &A, BigTorusVector &b);

void fixp_Ab_prod(BigTorusVector &res, BigTorusMatrix &A, BigTorusVector &b);

void fixp_Ab_prod_fake(BigTorusVector &res, BigTorusMatrix &A, BigTorusVector &b);

void fixp_tAb_prod_fake(BigTorusVector &res, BigTorusMatrix &A, BigTorusVector &b);

void fill_matrix_S(BigTorusMatrix &S);

void fill_matrix_Xy(BigTorusMatrix &X, BigTorusVector &y);

void fixp_sigmoid_vec(BigTorusVector &p, BigTorusVector &w, BigTorusVector &x);

void fixp_public_scale(BigTorusVector &res, int alpha);

void fixp_public_scale_fake(BigTorusVector &res, int alpha);

NTL::RR fixp_debug_norm(const BigTorusVector &v);

std::ostream &operator<<(std::ostream &out, const BigTorusVector &v);

uint8_t random_bit();

UINT64 random_uint64_t();

#endif //FHE_ARITHMETIC_H
