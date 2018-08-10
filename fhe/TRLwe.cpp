#include "TRLwe.h"

void pubKS(TRLwe &out, TLwe &in, pubKsKey &ks) {
    //const TRLweParams& out_params = out.params;
    const TLweParams& in_params = in.params;
    const int64_t out_prec_limbs = ks.out_prec_limbs;
    //const int64_t ks_prec_limbs = out_prec_limbs+2; //ks has 128-bit more precision
    const BigTorus& bitDecomp_in_offset = ks.bitDecomp_in_offset; // sum Bg/2 Bg^i
    __int128 bitDecomp_out_offset = ks.bitDecomp_out_offset; // -Bg/2

    BigTorus tmpDec(&in_params.fixp_params); //temp variable

    // out = trivial(b)
    trivial(out, in.getBT(), out_prec_limbs);
    for (UINT64 i=0; i<in_params.N; i++) {
        add(tmpDec, in.getAT(i), bitDecomp_in_offset);
        for (UINT64 j=0; j<ks.l_dec; ++j) {
            // coef aij of the decomposition
            __int128 aij = bitdecomp_coef128(tmpDec, j, out_prec_limbs) + bitDecomp_out_offset;
            // out = out - aij . ks_ij
            subMul(out, aij, ks.kskey[i][j], out_prec_limbs);
        }
    }
}
