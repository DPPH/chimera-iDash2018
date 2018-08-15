#include <cassert>
#include "BigComplex.h"

BigComplex::BigComplex(UINT64 nblimbs) : real(nblimbs), imag(nblimbs) {}

BigComplex::~BigComplex() {}

BigComplexRef::BigComplexRef(BigReal *real, BigReal *imag) : real(real), imag(imag) {}

BigComplexRef::BigComplexRef(const BigReal *real, const BigReal *imag) :
        real((BigReal *) real),
        imag((BigReal *) imag) {}

BigComplexRef::BigComplexRef(BigComplex &z) :
        real(&z.real),
        imag(&z.imag) {}

BigComplexRef::BigComplexRef(const BigComplex &z) :
        real((BigReal *) &z.real),
        imag((BigReal *) &z.imag) {}

void accurate_power_unity(BigComplexRef dest, int i, int n) {
    const UINT64 nblimbs = dest.real->nblimbs;
    const UINT64 FPREC = nblimbs * BITS_PER_LIMBS;
    mpfr_t angle;
    mpfr_t rsin;
    mpfr_t rcos;
    mpfr_init2(angle, FPREC + 1);
    mpfr_init2(rsin, FPREC + 1);
    mpfr_init2(rcos, FPREC + 1);
    //compute the angle
    i = ((i % n) + n) % n;
    mpfr_const_pi(angle, MPFR_RNDN);
    mpfr_mul_si(angle, angle, 2 * i, MPFR_RNDN);
    mpfr_div_si(angle, angle, n, MPFR_RNDN);
    //compute sin and cos
    mpfr_sin_cos(rsin, rcos, angle, MPFR_RNDN);
    //shift and write the answer
    mpfr_mul_2ui(rcos, rcos, FPREC, MPFR_RNDN);
    mpfr_mul_2ui(rsin, rsin, FPREC, MPFR_RNDN);
    mpfr_get_z(dest.real->value, rcos, MPFR_RNDN);
    mpfr_get_z(dest.imag->value, rsin, MPFR_RNDN);
    //cleanup
    mpfr_clear(angle);
    mpfr_clear(rsin);
    mpfr_clear(rcos);
}

void mul(BigComplexRef dest, BigComplexRef a, BigComplexRef b) {
    const UINT64 nblimbs = dest.real->nblimbs;
    BigReal t1(nblimbs);
    BigReal t2(nblimbs);
    BigReal t3(nblimbs);
    BigReal t4(nblimbs);
    mul(t1, *a.real, *b.real);
    mul(t2, *a.imag, *b.imag);
    mul(t3, *a.real, *b.imag);
    mul(t4, *a.imag, *b.real);
    sub(*dest.real, t1, t2);
    add(*dest.imag, t3, t4);
}

void mulTo(BigComplexRef dest, BigComplexRef a) {
    mul(dest, dest, a);
}

void add(BigComplexRef dest, BigComplexRef a, BigComplexRef b) {
    add(*dest.real, *a.real, *b.real);
    add(*dest.real, *a.real, *b.real);
    add(*dest.imag, *a.imag, *b.imag);
}

void sub(BigComplexRef dest, BigComplexRef a, BigComplexRef b) {
    sub(*dest.real, *a.real, *b.real);
    sub(*dest.imag, *a.imag, *b.imag);
}

void copy(BigComplexRef dest, BigComplexRef a) {
    copy(*dest.real, *a.real);
    copy(*dest.imag, *a.imag);
}

BigComplex *new_BigComplex_array(UINT64 n, UINT64 nblimbs) {
    BigComplex *reps = (BigComplex *) malloc(n * sizeof(BigComplex));
    for (UINT64 i = 0; i < n; i++) {
        new(reps + i) BigComplex(nblimbs);
    }
    return reps;
}

void delete_BigComplex_array(UINT64 n, BigComplex *array) {
    for (UINT64 i = 0; i < n; i++) {
        (array + i)->~BigComplex();
    }
    free(array);
}

void zero(BigComplexRef dest) {
    zero(*dest.real);
    zero(*dest.imag);
}

void addMulTo(BigComplexRef dest, const BigComplexRef &a, const BigComplexRef &b) {
    const UINT64 nblimbs = a.real->nblimbs;
    assert(dest.real->nblimbs == nblimbs);
    assert(dest.imag->nblimbs == nblimbs);
    assert(a.real->nblimbs == nblimbs);
    assert(a.imag->nblimbs == nblimbs);
    assert(b.real->nblimbs == nblimbs);
    assert(b.imag->nblimbs == nblimbs);
    BigComplex tmp(nblimbs);
    mul(tmp, a, b);
    add(dest, dest, tmp);
}

