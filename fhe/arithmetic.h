#ifndef FHE_ARITHMETIC_H
#define FHE_ARITHMETIC_H

#include "BigFixP.h"

void to_fixP(BigFixPRef reps, const NTL::RR &a);

void to_torus(BigTorusRef reps, const NTL::RR &a);

NTL::RR to_RR(const BigFixPRef &a);

void tAb_prod(BigFixPVector &res, BigFixPMatrix &A, BigFixPVector &b);

NTL::RR to_RR(const BigTorusRef &a);

void fill_matrix_S(BigFixPMatrix &S);

void fill_matrix_Xy(BigFixPMatrix &X, BigFixPVector &y);

void sigmoid_vec(BigFixPVector &p, BigFixPVector &w, BigFixPVector &x);

void public_scale(BigFixPVector &res, int alpha);

NTL::RR debug_norm(const BigFixPVector &v);

std::ostream &operator<<(std::ostream &out, const BigFixPVector &v);

#endif //FHE_ARITHMETIC_H
