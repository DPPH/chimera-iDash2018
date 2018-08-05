#include <complex>
#include <iostream>
#include <cstdint>
#include <vector>
#include <cassert>
#include <random>
#include <gmpxx.h>
#include <mpfr.h>

using namespace std;

typedef int64_t Torus64;

typedef int32_t Torus32;

#ifndef ALPHA_BITS
const int ALPHA_BITS = 60;
#endif

const int MIN_N = 1000 * ALPHA_BITS / 35;
const int LOG2N = int(ceil(log(MIN_N) / log(2.)));
const int N = 1 << LOG2N;
const int FPREC = ALPHA_BITS + 4;
const int IPREC = FPREC / 2;


struct Real96 {
    mpfr_t v;

    Real96() {
        mpfr_init2(v, FPREC + IPREC);
        mpfr_set_ui(v, 0, MPFR_RNDN);
    }

    Real96(const Real96 &r) {
        mpfr_init2(v, FPREC + IPREC);
        mpfr_set(v, r.v, MPFR_RNDN);
    }

    void operator=(const Real96 &r) { mpfr_set(v, r.v, MPFR_RNDN); }

    ~Real96() { mpfr_clear(v); }
};

void add(Real96 &dest, const Real96 &a, const Real96 &b) {
    mpfr_add(dest.v, a.v, b.v, MPFR_RNDN);
}

void sub(Real96 &dest, const Real96 &a, const Real96 &b) {
    mpfr_sub(dest.v, a.v, b.v, MPFR_RNDN);
}

void neg(Real96 &dest, const Real96 &a) {
    mpfr_neg(dest.v, a.v, MPFR_RNDN);
}

void extmul(Real96 &dest, int a, const Real96 &b) {
    mpfr_mul_si(dest.v, b.v, a, MPFR_RNDN);
}

void intmul(Real96 &dest, const Real96 &a, const Real96 &b) {
    mpfr_mul(dest.v, a.v, b.v, MPFR_RNDN);
}

void div2ui(Real96 &dest, const Real96 &a, uint64_t b) {
    mpfr_div_2ui(dest.v, a.v, b, MPFR_RNDN);
}

Real96 t64tor96(Torus64 v) {
    Real96 reps;
    int64_t vv = v;
    mpfr_set_si(reps.v, vv, MPFR_RNDN);
    mpfr_div_2ui(reps.v, reps.v, 64, MPFR_RNDN);
    return reps;
}

Real96 dtor96(double d) {
    Real96 reps;
    mpfr_set_d(reps.v, d, MPFR_RNDN);
    return reps;
}

//auto-deduced operators
Real96 operator+(const Real96 &a, const Real96 &b) {
    Real96 reps;
    add(reps, a, b);
    return reps;
}

void operator+=(Real96 &a, const Real96 &b) {
    add(a, a, b);
}

Real96 operator-(const Real96 &a, const Real96 &b) {
    Real96 reps;
    sub(reps, a, b);
    return reps;
}

void operator-=(Real96 &a, const Real96 &b) {
    sub(a, a, b);
}

Real96 operator-(const Real96 &a) {
    Real96 reps;
    neg(reps, a);
    return reps;
}

Real96 operator*(int a, const Real96 &b) {
    Real96 reps;
    extmul(reps, a, b);
    return reps;
}

void operator*=(Real96 &b, int a) {
    extmul(b, a, b);
}

Real96 operator*(const Real96 &a, const Real96 &b) {
    Real96 reps;
    intmul(reps, a, b);
    return reps;
}

void operator*=(Real96 &a, const Real96 &b) {
    intmul(a, a, b);
}


double r96tod(const Real96 &r) {
    return mpfr_get_d(r.v, MPFR_RNDN);
}

ostream &operator<<(ostream &out, const Real96 &r) {
    out << "[" << r96tod(r) << "]";
    return out;
}

bool very_close(const Real96 &a, const Real96 &b) {
    double dist = abs(r96tod(a - b));
    bool reps = dist < 1e-4;
    if (!reps) {
        cerr << "not close: " << a << " vs. " << b << endl;
    }
    return reps;
}





//-----------------------------------------------------------------------

using Cplx96=std::complex<Real96>;

//-----------------------------------------------------------------------


Real96 accurate_cos(int i, int n) { //cos(2pi*i/n)
    Real96 reps;
    mpfr_t res;
    mpfr_init2(res, FPREC + 1);
    i = ((i % n) + n) % n;
    mpfr_const_pi(res, MPFR_RNDN);
    mpfr_mul_si(res, res, 2 * i, MPFR_RNDN);
    mpfr_div_si(res, res, n, MPFR_RNDN);
    mpfr_cos(reps.v, res, MPFR_RNDN);
    return reps;
}

Real96 accurate_sin(int i, int n) { //sin(2pi*i/n)
    Real96 reps;
    mpfr_t res;
    mpfr_init2(res, FPREC + 1);
    i = ((i % n) + n) % n;
    mpfr_const_pi(res, MPFR_RNDN);
    mpfr_mul_si(res, res, 2 * i, MPFR_RNDN);
    mpfr_div_si(res, res, n, MPFR_RNDN);
    mpfr_sin(reps.v, res, MPFR_RNDN);
    return reps;
}

//reverse the bits of i (mod n)
int rev(int i, int n) {
    int reps = 0;
    for (int j = 1; j < n; j *= 2) {
        reps = 2 * reps + (i % 2);
        i /= 2;
    }
    return reps;
}

bool very_close(const Cplx96 &a, const Cplx96 &b) {
    bool reps = (very_close(a.real(), b.real()) && very_close(a.imag(), b.imag()));
    if (!reps) {
        cerr << "not close: " << a << " vs. " << b << endl;
    }
    return reps;
}


// FFT from Torus64^N to Cplx96^(N/2)  mod X^N+1
// N = 2048 (note: n=2N ici)


//at the beginning of iteration nn
// a_{j,i} has P_{i%nn}(omega^j) 
// where j between [rev(1) and rev(3)[
// and i between [0 and nn[
void ifft_check(int n, int nn, const Cplx96 *acur, const vector<Real96> &a, const vector<Cplx96> &powomega) {
    int ns4 = n / 4;
    cerr << "Checking iteration " << nn << endl;
    for (int i = 0; i < ns4; i++) {
        cout << "i: " << i << "   " << acur[i] << endl;
    }
    int m = n / nn;
    int rev1m = rev(1, m);
    int rev3m = rev(3, m);
    int idex = 0;
    for (int revj = rev1m; revj < rev3m; revj++) {
        int j = rev(revj, m);
        cerr << "check-- j: " << j << endl;
        for (int i = 0; i < nn; i++) {
            cerr << "check--- i: " << i << "(mod " << nn << ")" << endl;
            const Cplx96 &test_cur = acur[idex];
            //sum_[t=i%nn] a_t omega^jt
            Cplx96 pij(t64tor96(0), t64tor96(0));
            for (int k = i; k < n; k += nn) {
                cout << "a_" << k << ": " << a[k] << endl;
                pij += a[k] * powomega[(k * j) % n];
            }
            assert(very_close(test_cur, pij));
            idex++;
        }
    }
}


//at the beginning of iteration halfnn:
//   m=n/halfnn
//   P_{j%m}(omb^i)
//   for j in [rev(1,m) to rev(3,m)[
//   and i in [0,halfnn[
void fft_check(
        int n, int halfnn,
        const Cplx96 *pcur,
        const vector<Cplx96> &p,
        const vector<Cplx96> &powombar
) {
    int ns4 = n / 4;
    cerr << "DIRECT FFT: Checking iteration " << halfnn << endl;
    for (int i = 0; i < ns4; i++) {
        cout << "i: " << i << "   " << pcur[i] << endl;
    }
    int m = n / halfnn;
    int rev1m = rev(1, m);
    int rev3m = rev(3, m);
    int idex = 0;
    for (int revj = rev1m; revj < rev3m; revj++) {
        int j = rev(revj, m);
        cerr << "check-- j: " << j << "(mod " << m << ")" << endl;
        for (int i = 0; i < halfnn; i++) {
            cerr << "check--- i: " << i << endl;
            //P_sum_[k=j%m] p_k omb^ik-j
            Cplx96 pij(t64tor96(0), t64tor96(0));
            for (int k = j; k < n; k += m) {
                //if (halfnn==8 && j==1 && i==1) cerr << "pij(" << pij << ")" << "+= p_"<<k<<"("<<p[k]<<") * omb["<<i*(k-j)<<"]("<< powombar[(i*(k-j)) % n] <<")" << endl;
                pij += p[k] * powombar[(i * (k - j)) % n];
            }
            assert(very_close(pcur[idex], pij));
            idex++;
        }
    }
}


void precomp_iFFT(vector<Cplx96> &powomega, int n) {
    powomega.resize(n);
    for (int i = 0; i < n; i++)
        powomega[i] = Cplx96(accurate_cos(i, n), accurate_sin(i, n));
}

void precomp_FFT(vector<Cplx96> &powombar, int n) {
    powombar.resize(n);
    for (int i = 0; i < n; i++)
        powombar[i] = Cplx96(accurate_cos(i, n), accurate_sin((n - i) % n, n));
}

// P -> P(omega)
void iFFT(Cplx96 *out, const Real96 *in, int n, const vector<Cplx96> &powomega) {
    //const int N = n/2;
    const int ns4 = n / 4;

#ifndef NDEBUG
    vector<Real96> a; a.resize(n);
    for (int i=0; i<n/2; i++)
        div2ui(a[i],in[i],1);
    for (int i=0; i<n/2; i++)
        a[n/2+i]=-a[i];
#endif


    //interpret the input coefs as real and imaginary parts:
    //RRRRRRRRRRIIIIIIIII
    //multiply by omega^j
    for (int j = 0; j < ns4; j++)
        out[j] = Cplx96(in[j], in[j + ns4]) * powomega[j];

    //at the beginning of iteration nn
    // a_{j,i} has P_{i%nn}(omega^j) 
    // where j between [rev(1) and rev(3)[
    // and i between [0 and nn[
    for (int nn = ns4; nn >= 2; nn /= 2) {
        int halfnn = nn / 2;
#ifndef NDEBUG
        cerr << "Starting iteration " << nn << endl;
        int m = n/nn;
        ifft_check(n, nn, out, a, powomega);
#endif
        for (int block = 0; block < ns4; block += nn) {
#ifndef NDEBUG
            int j = rev(rev(1,m)+block,m);
            cerr << "-- block j: " << j  << " --> " << j << "," << j+halfnn/2 << endl;
#endif
            for (int off = 0; off < halfnn; off++) {
#ifndef NDEBUG
                cerr << "--- i: " << off << " using: omg^" << (2*(ns4/halfnn)*off)%n << endl;
#endif
                Cplx96 t1 = out[block + off];
                Cplx96 t2 = out[block + off + halfnn];
                out[block + off] = t1 + t2;
                out[block + off + halfnn] = (t1 - t2) * powomega[(2 * (ns4 / halfnn) * off) % n];
            }
        }
    }
    {
#ifndef NDEBUG
        int nn = 1;
        ifft_check(n, nn, out, a, powomega);
#endif
    }
}

// P(omega) -> P
void FFT(Real96 *out, Cplx96 *in, int n, const vector<Cplx96> &powombar) {
    //const int N = n/2;
    const int ns4 = n / 4;

#ifndef NDEBUG
    vector<Cplx96> a; a.resize(n);
    for (int i=0; i<n; i++) a[i]=Cplx96(t64tor96(0),t64tor96(0));
    int rev1m = rev(1,n);
    int rev3m = rev(3,n);
    for (int revj=rev1m; revj<rev3m; revj++) {
        int j = rev(revj,n);
        cout << "assign:" << j << " " << revj-rev1m << endl;
        a[j]=in[revj-rev1m];
        a[n-j]=Cplx96(a[j].real(),-a[j].imag());
    }
#endif


#ifndef NDEBUG
    cerr << "Checking iteration 1" << endl;
    fft_check(n, 1, in, a, powombar);
#endif

    //at the beginning of iteration nn
    // a_{j,i} has P_{i%nn}(omega^j) 
    // where j between [rev(1) and rev(3)[
    // and i between [0 and nn[
    for (int nn = 2; nn <= ns4; nn *= 2) {
        int halfnn = nn / 2;
        for (int block = 0; block < ns4; block += nn) {
            //#ifndef NDEBUG
            //	    int j = rev(rev(1,m)+block,m);
            //	    cerr << "-- block j: " << j  << " --> " << j << "," << j+halfnn/2 << endl;
            //#endif
            for (int off = 0; off < halfnn; off++) {
                //#ifndef NDEBUG
                //		cerr << "--- i: " << off << " using: omg^" << (2*(ns4/halfnn)*off)%n << endl;
                //#endif
                Cplx96 t1 = in[block + off];
                Cplx96 t2 = in[block + off + halfnn] * powombar[(2 * (ns4 / halfnn) * off) % n];
                in[block + off] = t1 + t2;
                in[block + off + halfnn] = (t1 - t2);
            }
        }
#ifndef NDEBUG
        cerr << "Ending iteration " << nn << endl;
        fft_check(n, nn, in, a, powombar);
#endif
    }


    //interpret the input coefs as real and imaginary parts:
    //RRRRRRRRRRIIIIIIIII
    //multiply by omega^j
    for (int j = 0; j < ns4; j++) {
        in[j] *= powombar[j];
        div2ui(out[j], in[j].real(), LOG2N - 1);      // /ns4;  //divide by N/2
        div2ui(out[j + ns4], in[j].imag(), LOG2N - 1);  // /ns4; //divide by N/2
    }

    {
        //#ifndef NDEBUG
        //	int nn = 1;
        //	ifft_check(n, nn, out, a, powomega);
        //#endif
    }
}


int main(int argc, char **argv) {

    const int Ns2 = N / 2;
#ifndef NDEBUG
    const int NBTRIALS=1;
#else
    const int NBTRIALS = 100;
#endif

#ifndef NDEBUG
    cout << "trigo_test" << endl;
    //test cosinus
    for (int i=0; i<=2*N; i++) {
        //if (i!=1024) continue;
        Real96 c = accurate_cos(i,2*N);
        Real96 s = accurate_sin(i,2*N);
        Real96 t = c*c+s*s;
        if (!very_close(t, dtor96(1))) {
            cout << "i: " << i  << ", " << cos(2*i*M_PI/64.) << endl;
            cout << "c " << c << endl;
            cout << "s " << s << endl;
            cout << "c*c " << c*c << endl;
            cout << "s*s " << s*s << endl;
            cout << "cos test" <<  i << ":" << t << endl;
            assert(false );
        }
    }
    cout << "...passed!" << endl;
#endif

    //test de la IFFT
    Cplx96 out[Ns2];
    Real96 in[N];
    Real96 revout[N];

    vector<Cplx96> powomega;
    vector<Cplx96> powombar;

    std::default_random_engine generator;
    std::uniform_int_distribution<Torus64> distribution(numeric_limits<Torus64>::min(), numeric_limits<Torus64>::max());

    for (int i = 0; i < N; i++) {
        in[i] = t64tor96(distribution(generator));
    }

    precomp_iFFT(powomega, 2 * N);
    precomp_FFT(powombar, 2 * N);

#ifndef NDEBUG
    //test powomega
    cout << "powomega_test" << endl;
    for (int i=0; i<2*N; i++) {
        //if (i!=1024) continue;
        //cout << "powomega: " << i << " : " << powomega[i] << endl;
        assert(very_close(powomega[i]*powomega[(2*N-i)%(2*N)],Cplx96(dtor96(1))));
        assert(very_close(powomega[i]*powombar[i],Cplx96(dtor96(1))));
    }
    cout << "...passed!" << endl;
#endif

    clock_t t0 = clock();
    for (int i = 0; i < NBTRIALS; i++)
        iFFT(out, in, 2 * N, powomega);
    clock_t t1 = clock();
    for (int i = 0; i < NBTRIALS; i++)
        FFT(revout, out, 2 * N, powombar);
    clock_t t2 = clock();
#ifndef NDEBUG
    for (int i=0; i<N; i++)
        cout << hex << revout[i] << " "<< in[i] << " " << (revout[i]-in[i]) << endl;
#endif
    cerr << "ALPHA_BITS: " << ALPHA_BITS << endl;
    cerr << "MIN_N: " << MIN_N << endl;
    cerr << "LOG2N: " << LOG2N << endl;
    cerr << "N: " << N << endl;
    cerr << "FPREC: " << FPREC << endl;
    cerr << "IPREC: " << IPREC << endl;
    cerr << "time_IFFT: " << (t1 - t0) / double(NBTRIALS) << " mus" << endl;
    cerr << "time_FFT: " << (t2 - t1) / double(NBTRIALS) << " mus" << endl;
    printf("%d %d %d %d %d %d %lf %lf\n", ALPHA_BITS, MIN_N, LOG2N, N, FPREC, IPREC, (t1 - t0) / double(NBTRIALS),
           (t2 - t1) / double(NBTRIALS));

}
