#include <benchmark/benchmark.h>
#include "../BigReal.h"
#include "../BigComplex.h"
#include "../BigFFT.h"

NTL_CLIENT

static void RangesAll(benchmark::internal::Benchmark* b) {
    for (int n = 9; n < 16; ++n) {
        b->Args({1<<n, 1});
        for (int nblimbs = 5; nblimbs < 20; nblimbs += 5)
          b->Args({1<<n, nblimbs});
    }
}

class FFT_Bench : public benchmark::Fixture {
public:
    void SetUp(const ::benchmark::State &state) {
        n = state.range(0);
        const auto nblimbs = state.range(1);

        N = n / 2;
        ns4 = n / 4;

        powombar = precomp_FFT(n, nblimbs);

        inp = new_BigComplex_array(ns4, nblimbs);
        out = new_BigReal_array(N, nblimbs);

        int64_t PREC = nblimbs * BITS_PER_LIMBS;
        RR::SetPrecision(PREC);
        vec_RR test_in;
        test_in.SetLength(N);
        // fill random values
        for (int i = 0; i < ns4; i++) {
            test_in[2 * i] = random_RR() / ns4 / ns4;
            test_in[2 * i + 1] = random_RR() / ns4 / ns4;
            to_BigReal(inp[i].real, test_in[2 * i]);
            to_BigReal(inp[i].imag, test_in[2 * i + 1]);
        }
    }

    void TearDown(const ::benchmark::State &state) {
        //clean up
        delete_BigReal_array(N, out);
        delete_BigComplex_array(ns4, inp);

        clear_precomp_FFT(powombar);
    }

    int n;
    int N;
    int ns4;

    BigComplex *powombar;
    BigComplex *inp;
    BigReal *out;
};

BENCHMARK_DEFINE_F(FFT_Bench, FFT)(benchmark::State &state) {
    int64_t a = rand();
    for (auto _ : state) {
        FFT(out, inp, n, powombar);
    }
}
BENCHMARK_REGISTER_F(FFT_Bench, FFT)
    // ->Ranges({  {1<<9, 1<<16},
    //             {1, 16} })
    ->Apply(RangesAll)
    ->Unit(benchmark::kMicrosecond);



class iFFT_Bench : public benchmark::Fixture {
public:
    void SetUp(const ::benchmark::State &state) {
        n = state.range(0);
        const auto nblimbs = state.range(1);

        N = n / 2;
        ns4 = n / 4;

        powomega = precomp_iFFT(n, nblimbs);

        inp = new_BigReal_array(N, nblimbs);
        out = new_BigComplex_array(ns4, nblimbs);

        int64_t PREC = nblimbs * BITS_PER_LIMBS;
        RR::SetPrecision(PREC);
        vec_RR test_in;
        test_in.SetLength(N);
        // fill random values
        for (int i = 0; i < N; i++) {
            test_in[i] = random_RR() / ns4 / ns4;
            to_BigReal(inp[i], test_in[i]);
        }
    }

    void TearDown(const ::benchmark::State &state) {
        //clean up
        delete_BigReal_array(N, inp);
        delete_BigComplex_array(ns4, out);

        clear_precomp_iFFT(powomega);
    }

    int n;
    int N;
    int ns4;

    BigComplex *powomega;
    BigReal *inp;
    BigComplex *out;
};



BENCHMARK_DEFINE_F(iFFT_Bench, iFFT)(benchmark::State &state) {
    int64_t a = rand();
    for (auto _ : state) {
        iFFT(out, inp, n, powomega);
    }
}
BENCHMARK_REGISTER_F(iFFT_Bench, iFFT)
    // ->Ranges({  {1<<9, 1<<16},
    //             {1, 16} })
    ->Apply(RangesAll)
    ->Unit(benchmark::kMicrosecond);
