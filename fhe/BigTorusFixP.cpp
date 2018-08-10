#include <gmp.h>
#include <cassert>
#include "BigTorus.h"

fixp_add_params::fixp_add_params() :
        ia_limbs(0), ib_limbs(0), ic_limbs(0),
        oc_limbs(0), sa(0), sb(0),
        tmp_a(nullptr), tmp_b(nullptr) {}

void fixp_prepareAdd(fixp_add_params &out, const BigTorusParams &pres, const BigTorusParams &pa,
                     const BigTorusParams &pb,
                     int64_t out_precision_bits) {
    out.ia_limbs = pa.torus_limbs;
    out.ib_limbs = pb.torus_limbs;
    out.ic_limbs = pres.torus_limbs;
    out.oc_limbs = (out_precision_bits + BITS_PER_LIMBS - 1) / BITS_PER_LIMBS;
    assert_dramatically(out.oc_limbs <= out.ic_limbs,
                        "too much precision required"
    );

    // verify (weak) rho_c >= max(rho_a, rho_b)
    int64_t minimal_plaintext_expo = std::max(pa.plaintext_expo, pb.plaintext_expo);
    assert_weakly(pres.plaintext_expo >= minimal_plaintext_expo,
                  "addition plaintext exponent is too small"
    );

    // verify (strong) L_c <= min(rho_a + L_a, rho_b + L_b) - rho_c
    int64_t maximal_level_expo = std::min(pa.plaintext_expo + pa.level_expo, pb.plaintext_expo + pb.level_expo)
                                 - pres.plaintext_expo;
    assert_dramatically(pres.level_expo <= maximal_level_expo,
                        "addition level exponent is too small"
    );

    // left shift computation for a
    out.sa = pa.plaintext_expo + pa.level_expo - pres.plaintext_expo - pres.level_expo;
    // left shift computation for b
    out.sb = pb.plaintext_expo + pb.level_expo - pres.plaintext_expo - pres.level_expo;

    //allocate tmp buffer
    out.tmp_a = new UINT64[out.oc_limbs + 1];
    out.tmp_b = new UINT64[out.oc_limbs + 1];
}

void fixp_releaseAdd(fixp_add_params &out) {
    delete[] out.tmp_a;
    delete[] out.tmp_b;
}

void fixp_special_add_lshift(UINT64 *out, UINT64 *a, int64_t out_limbs, int64_t a_limbs, int64_t bit_shift) {
    int64_t sbits = bit_shift % BITS_PER_LIMBS;
    int64_t slimbs = bit_shift / BITS_PER_LIMBS;
    int64_t oc = out_limbs;
    int64_t signif_limbs = a_limbs - slimbs;

    if (sbits == 0) {
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
        int64_t ocp = oc + 1;
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

void fixp_raw_add(UINT64 *reps, UINT64 *a, UINT64 *b, const fixp_add_params &params) {
    UINT64 *tmp_a = params.tmp_a;
    UINT64 *tmp_b = params.tmp_b;
    UINT64 *c = reps + params.ic_limbs - params.oc_limbs;

    //left shift a by sa positions and keep oc limbs
    fixp_special_add_lshift(tmp_a, a, params.oc_limbs, params.ia_limbs, params.sa);
    //left shift b by sb positions and keep oc limbs
    fixp_special_add_lshift(tmp_b, b, params.oc_limbs, params.ib_limbs, params.sb);
    //add a,b (on oc_limbs)
    mpn_add_n(c, tmp_a, tmp_b, params.oc_limbs);
}

void fixp_raw_sub(UINT64 *reps, UINT64 *a, UINT64 *b, const fixp_add_params &params) {
    UINT64 *tmp_a = params.tmp_a;
    UINT64 *tmp_b = params.tmp_b;
    UINT64 *c = reps + params.ic_limbs - params.oc_limbs;

    //left shift a by sa positions and keep oc limbs
    fixp_special_add_lshift(tmp_a, a, params.oc_limbs, params.ia_limbs, params.sa);
    //left shift b by sb positions and keep oc limbs
    fixp_special_add_lshift(tmp_b, b, params.oc_limbs, params.ib_limbs, params.sb);
    //add a,b (on oc_limbs)
    mpn_sub_n(c, tmp_a, tmp_b, params.oc_limbs);
}


void fixp_add(BigTorusRef reps, const BigTorusRef &a, const BigTorusRef &b, UINT64 out_precision_bits) {
    fixp_add_params addParams;
    if (out_precision_bits == NA) out_precision_bits = reps.params->torus_limbs * BITS_PER_LIMBS;
    fixp_prepareAdd(addParams, *reps.params, *a.params, *b.params, out_precision_bits);
    fixp_raw_add(reps.limbs, a.limbs, b.limbs, addParams);
    fixp_releaseAdd(addParams);
}

void fixp_sub(BigTorusRef reps, const BigTorusRef &a, const BigTorusRef &b, UINT64 out_precision_bits) {
    fixp_add_params addParams;
    if (out_precision_bits == NA) out_precision_bits = reps.params->torus_limbs * BITS_PER_LIMBS;
    fixp_prepareAdd(addParams, *reps.params, *a.params, *b.params, out_precision_bits);
    fixp_raw_sub(reps.limbs, a.limbs, b.limbs, addParams);
    fixp_releaseAdd(addParams);
}

