#ifndef FHE_ARITHMETIC_H
#define FHE_ARITHMETIC_H

#include "BigFixP.h"
#include "BigFixPVector.h"

void to_fixP(BigFixPRef reps, const NTL::RR &a);

NTL::RR to_RR(const BigFixPRef &a);

void tAb_prod(BigFixPVector &res, BigFixPMatrix &A, BigFixPVector &b);

void Ab_prod(BigFixPVector &res, BigFixPMatrix &A, BigFixPVector &b);

void Ab_prod_fake(BigFixPVector &res, BigFixPMatrix &A, BigFixPVector &b);

void tAb_prod_fake(BigFixPVector &res, BigFixPMatrix &A, BigFixPVector &b);

void fill_matrix_S(BigFixPMatrix &S);

void fill_matrix_Xy(BigFixPMatrix &X, BigFixPVector &y);

void sigmoid_vec(BigFixPVector &p, BigFixPVector &w, BigFixPVector &x);

void public_scale(BigFixPVector &res, int alpha);

void public_scale_fake(BigFixPVector &res, int alpha);

NTL::RR debug_norm(const BigFixPVector &v);

std::ostream &operator<<(std::ostream &out, const BigFixPVector &v);

uint8_t random_bit();

uint64_t random_uint64_t();

#endif //FHE_ARITHMETIC_H
