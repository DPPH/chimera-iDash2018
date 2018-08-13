#include <benchmark/benchmark.h>
#include <commons.h>
#include <BigTorus.h>
#include <TLwe.h>
#include <TRLwe.h>

using namespace std;

static void RangesAll(benchmark::internal::Benchmark *b) {
    for (int n = 32; n <= 16384; n *= 2) {
        b->Args({n});
    }
}

std::shared_ptr<pubKsKey>
ks_keygen_fake(const TRLweParams &out_params, const TLweParams &in_params,
               const TLweKey &in_key, const TLweKey &out_key,
               const UINT64 out_alpha_bits) {
    //automatic parameter deduction
    const double out_variance_bits = out_alpha_bits * 2;
    const double in_variance_bits = out_variance_bits + 1;     // 1/2 from the input noise
    const double errdec_variance_bits = out_variance_bits + 2; // 1/4 from decomp error
    const double kssum_variance_bits = out_variance_bits + 2;  // 1/4 from the big sum


    //smallest decomposition in base 2^128 giving an error variance <= errdec_variance_bits
    const UINT64 ldec = ceil(errdec_variance_bits / 2. / 128.);

    //the input noise variance must be smaller than in_variance_bits
    // -> input nblimbs must be larger than limb_prec(variance/2)
    assert_dramatically(in_params.fixp_params.torus_limbs >= UINT64(limb_precision(in_variance_bits / 2)));

    //smallest variance for the keyswitch key
    const double ks_variance_bits = kssum_variance_bits + log2(in_params.N) + log2(ldec) + 2 * 127.;
    const UINT64 ks_alpha_bits = UINT64(ks_variance_bits / 2);
    const UINT64 ks_nblimbs = UINT64(limb_precision(ks_alpha_bits) + 1) & UINT64(-2); //must be even
    TRLweParams ks_params(out_params.N, ks_nblimbs);

    //the output nblimbs must be larger than limb_prec(out_alpha_bits)
    assert_dramatically(out_params.fixp_params.torus_limbs >= UINT64(limb_precision(out_alpha_bits)));

    pubKsKey *reps = new pubKsKey(in_params, out_params, TRLweParams(ks_params), ldec);

    //plaintext must have the same precision as the keyswitch key
    BigTorusPolynomial plaintext(ks_params.N, ks_params.fixp_params);
    for (UINT64 i = 0; i < in_params.N; i++) {
        for (UINT64 j = 1; j <= ldec; j++) {
            // plaintext = ks[i] / Bg^j
            zero(plaintext);
            if (in_key.key[i] == 1) { plaintext.getAT(0).limbs_end[-2 * j] = 1; }
            //native_encrypt(reps->kskey[i][j - 1], plaintext, out_key, ks_alpha_bits);
            random(reps->kskey[i][j - 1].a[0], ks_params.fixp_params.torus_limbs);
            random(reps->kskey[i][j - 1].a[1], ks_params.fixp_params.torus_limbs);
        }
    }

    // sum Bg/2 Bg^-i
    zero(reps->bitDecomp_in_offset);
    for (UINT64 i = 1; i <= reps->in_params.fixp_params.torus_limbs; i += 2) {
        reps->bitDecomp_in_offset.limbs_end[-i] = 0x8000000000000000UL;
    }

    reps->bitDecomp_out_offset = (__int128(1) << __int128(127)); // -Bg/2
    return std::shared_ptr<pubKsKey>(reps);
}

class pubKs_Bench : public benchmark::Fixture {
public:

    unique_ptr<BigTorusParams> bt_params_in;
    unique_ptr<BigTorusParams> bt_params_out;
    unique_ptr<TLweParams> tlwe_params_in;
    unique_ptr<TRLweParams> trlwe_params_out;
    shared_ptr<TLweKey> key_in;
    shared_ptr<TLweKey> key_out;
    shared_ptr<pubKsKey> ks_key;
    int64_t limb_prec;
    unique_ptr<TLwe> ciphertext;
    unique_ptr<TRLwe> res;

    void SetUp(const ::benchmark::State &state) {
        int64_t N_in = state.range(0);
        int64_t N_out = 2 * N_in;
        int64_t nblimbs_in = 5;
        int64_t nblimbs_out = 6;
        int64_t alpha_bits = nblimbs_in * BITS_PER_LIMBS - 2;

        limb_prec = limb_precision(alpha_bits);
        bt_params_in.reset(new BigTorusParams(nblimbs_in));
        bt_params_out.reset(new BigTorusParams(nblimbs_out));
        tlwe_params_in.reset(new TLweParams(N_in, *bt_params_in));
        trlwe_params_out.reset(new TRLweParams(N_out, *bt_params_out));

        key_in = tlwe_keygen(*tlwe_params_in);
        key_out = tlwe_keygen(*trlwe_params_out);


        ks_key = ks_keygen_fake(
                *trlwe_params_out, *tlwe_params_in,
                *key_in, *key_out, alpha_bits);


        // take a random plaintext
        BigTorus plaintext(*bt_params_in);
        random(plaintext);

        // encrypt it
        ciphertext.reset(new TLwe(*tlwe_params_in));
        native_encrypt(*ciphertext, plaintext, *key_in);

        res.reset(new TRLwe(*trlwe_params_out));

    }

    void TearDown(const ::benchmark::State &state) {

        ks_key = nullptr;
        ciphertext = nullptr;
        res = nullptr;
    }


};


BENCHMARK_DEFINE_F(pubKs_Bench, pubKS)(benchmark::State &state) {
    for (auto _ : state) {
        pubKS(*res, *ciphertext, *ks_key, limb_prec);
    }
}


BENCHMARK_REGISTER_F(pubKs_Bench, pubKS)

        ->Apply(RangesAll)
        ->Unit(benchmark::kMicrosecond);
