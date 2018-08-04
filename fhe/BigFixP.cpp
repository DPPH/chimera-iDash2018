#include "BigFixP.h"

BigFixPAddParams::BigFixPAddParams() :
        ia_limbs(0), ib_limbs(0), ic_limbs(0),
        oc_limbs(0), sa(0), sb(0),
        tmp_a(nullptr), tmp_b(nullptr) {}

void prepareAdd(BigFixPAddParams &out, const BigFixPParams &pres, const BigFixPParams &pa, const BigFixPParams &pb,
                int64_t out_precision_bits) {
    out.ia_limbs = pa.torus_params.torus_limbs;
    out.ib_limbs = pb.torus_params.torus_limbs;
    out.ic_limbs = pres.torus_params.torus_limbs;
    out.oc_limbs = (out_precision_bits + BITS_PER_LIMBS - 1) / BITS_PER_LIMBS;
    assert_dramatically(out.oc_limbs <= out.ic_limbs);

    // verify (weak) rho_c >= max(rho_a, rho_b)
    assert_weakly(pres.plaintext_expo >= std::max(pa.plaintext_expo, pb.plaintext_expo));

    // verify (strong) L_c <= min(rho_a + L_a, rho_b + L_b) - rho_c
    assert_dramatically(pres.level_expo <=
                        std::min(pa.plaintext_expo + pa.level_expo, pb.plaintext_expo + pb.level_expo)
                        - pres.plaintext_expo);

    // left shift computation for a
    out.sa = pa.plaintext_expo + pa.level_expo - pres.plaintext_expo - pres.level_expo;
    // left shift computation for b
    out.sb = pb.plaintext_expo + pb.level_expo - pres.plaintext_expo - pres.level_expo;

    //allocate tmp buffer
    out.tmp_a = new uint64_t[out.oc_limbs + 1];
    out.tmp_b = new uint64_t[out.oc_limbs + 1];
}

void releaseAdd(BigFixPAddParams &out) {
    delete[] out.tmp_a;
    delete[] out.tmp_b;
}

void special_add_lshift(uint64_t *out, uint64_t *a, int64_t out_limbs, int64_t a_limbs, int64_t bit_shift) {
    int64_t sbits = bit_shift % BITS_PER_LIMBS;
    int64_t slimbs = bit_shift / BITS_PER_LIMBS;
    int64_t oc = out_limbs;
    int64_t signif_limbs = a_limbs - slimbs;

    if (sbits==0) {
        //it is a shift by an exact amount of limbs
        if (signif_limbs <= 0) {
            // a:        **************|
            // +slimbs:  ----------------->
            //out is all zero
            mpn_zero(out, oc);
        } else {
            assert(signif_limbs > 0);
            if (signif_limbs >= oc) {
                // a:        **************|
                // +slimbs:  ----->********XXXXXX|
                // oc:              _______|
                //copy oc values from pos signif-oc to signif
                mpn_copyi(out, a + signif_limbs - oc, oc);
            } else {
                // a:        **************|
                // +slimbs:  ----->********XXXXXX|
                // oc:         ____________|
                // put oc-(signif) zeros on the left
                // copy signif values on the right
                mpn_zero(out, oc - signif_limbs);
                mpn_copyi(out + oc - signif_limbs, a, signif_limbs);
            }
        }
    } else {
        int64_t ocp = oc+1;
        //we shift a by slimbs limbs + sbits bits
        //and we keep oc+1 limbs of the output
        //right shift by 64-sbits bits if not exact
        //it is a shift by an exact amount of limbs
        if (signif_limbs <= 0) {
            // a:        **************|
            // +slimbs:  ----------------->
            //out is all zero
            mpn_zero(out, oc);
        } else {
            assert(signif_limbs > 0);
            if (signif_limbs >= ocp) {
                // a:        **************|
                // +slimbs:  ----->********XXXXXX|
                // oc:             -_______|
                //copy oc values from pos signif-oc to signif
                mpn_copyi(out, a + signif_limbs - ocp, ocp);
            } else {
                // a:        **************|
                // +slimbs:  ----->********XXXXXX|
                // oc:        -____________|
                // put ocp-(signif) zeros on the left
                // copy signif values on the right
                mpn_zero(out, ocp - signif_limbs);
                mpn_copyi(out + ocp - signif_limbs, a, signif_limbs);
            }
            //finally, right shift by 64-sbits
            mpn_rshift(out, out, ocp, BITS_PER_LIMBS - sbits);
            out[out_limbs] = 0;
        }
    }
}

void fixPRawAdd(uint64_t *reps, uint64_t *a, uint64_t *b, const BigFixPAddParams &params) {
    uint64_t *tmp_a = params.tmp_a;
    uint64_t *tmp_b = params.tmp_b;
    uint64_t *c = reps + params.ic_limbs - params.oc_limbs;

    //left shift a by sa positions and keep oc limbs
    special_add_lshift(tmp_a, a, params.oc_limbs, params.ia_limbs, params.sa);
    //left shift b by sb positions and keep oc limbs
    special_add_lshift(tmp_b, b, params.oc_limbs, params.ib_limbs, params.sb);
    //add a,b (on oc_limbs)
    mpn_add_n(c, tmp_a, tmp_b, params.oc_limbs);
}

void add(BigFixP &reps, const BigFixP &a, const BigFixP &b, int64_t out_precision_bits) {
    BigFixPAddParams addParams;
    if (out_precision_bits == NA) out_precision_bits = reps.params->torus_params.torus_limbs * BITS_PER_LIMBS;
    prepareAdd(addParams, *reps.params, *a.params, *b.params, out_precision_bits);
    fixPRawAdd(reps.limbs_raw, a.limbs_raw, b.limbs_raw, addParams);
    releaseAdd(addParams);
}

BigTorusRef::BigTorusRef(uint64_t *limbs_raw, const BigTorusParams *params) : limbs_raw(limbs_raw), params(params) {}

BigTorusRef::BigTorusRef(const BigTorus &torus): limbs_raw(torus.limbs_raw), params(torus.params) {}

BigTorusRef::BigTorusRef(BigTorus &torus): limbs_raw(torus.limbs_raw), params(torus.params) {}

BigFixPRef::BigFixPRef(const BigFixP &a): limbs_raw(a.limbs_raw), params(a.params) {}
