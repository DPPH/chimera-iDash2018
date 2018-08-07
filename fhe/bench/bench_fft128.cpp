#include <benchmark/benchmark.h>
#include "../BigFFT128.h"

static void RangesAll(benchmark::internal::Benchmark *b) {
    for (int n = 9; n < 16; ++n) {
        b->Args({1 << n, 1});
    }
}

class FFT128_Bench : public benchmark::Fixture {
public:
    void SetUp(const ::benchmark::State &state) {
        n = state.range(0);
        const auto nblimbs = state.range(1);

        N = n / 2;
        ns4 = n / 4;

        powombar = precomp_FFT128(n, nblimbs);

        inp = new Complex128[ns4];
        out = new Real128[N];

        // fill random values
        for (int i = 0; i < ns4; i++) {
            inp[i] = Complex128(rand() / double(RAND_MAX), rand() / double(RAND_MAX));
            //inp[i] = Complex128 (to_RR(rand() / double(RAND_MAX)), to_RR(rand() / double(RAND_MAX)));
        }
    }

    void TearDown(const ::benchmark::State &state) {
        //clean up
        delete[] out;
        delete[] inp;

        clear_precomp_FFT128(powombar);
    }

    int n;
    int N;
    int ns4;

    Complex128 *powombar;
    Complex128 *inp;
    Real128 *out;
};

BENCHMARK_DEFINE_F(FFT128_Bench, FFT
)(
benchmark::State &state
) {
for (
auto _
: state) {
FFT(out, inp, n, powombar
);
}
}
BENCHMARK_REGISTER_F(FFT128_Bench, FFT
)
// ->Ranges({  {1<<9, 1<<16},
//             {1, 16} })
->
Apply(RangesAll)
->
Unit(benchmark::kMicrosecond);


class iFFT128_Bench : public benchmark::Fixture {
public:
    void SetUp(const ::benchmark::State &state) {
        n = state.range(0);
        const auto nblimbs = state.range(1);

        N = n / 2;
        ns4 = n / 4;

        powomega = precomp_iFFT128(n, nblimbs);

        inp = new Real128[N];
        out = new Complex128[ns4];

        // fill random values
        for (int i = 0; i < N; i++) {
            inp[i] = rand() / double(RAND_MAX);
        }
    }

    void TearDown(const ::benchmark::State &state) {
        //clean up
        delete[] inp;
        delete[] out;

        clear_precomp_iFFT128(powomega);
    }

    int n;
    int N;
    int ns4;

    Complex128 *powomega;
    Real128 *inp;
    Complex128 *out;
};


BENCHMARK_DEFINE_F(iFFT128_Bench, iFFT
)(
benchmark::State &state
) {
for (
auto _
: state) {
iFFT(out, inp, n, powomega
);
}
}
BENCHMARK_REGISTER_F(iFFT128_Bench, iFFT
)
// ->Ranges({  {1<<9, 1<<16},
//             {1, 16} })
->
Apply(RangesAll)
->
Unit(benchmark::kMicrosecond);
