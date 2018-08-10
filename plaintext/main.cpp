#include <io.h>
#include <common.h>

#include <NTL/LLL.h>
#include <NTL/ZZ_limbs.h>
#include <cassert>

NTL_CLIENT;

typedef NTL::Vec<float> vec_float;
typedef NTL::Mat<float> mat_float;


vec_float sigmoid_vec(const vec_float &x) {
    const int n = x.length();
    vec_float reps;
    reps.SetLength(n);
    for (int i = 0; i < n; i++) reps[i] = sigmoid(x[i]);
    return reps;
}

void clear(vec_float &v) {
    const int n = v.length();
    for (int i = 0; i < n; i++) v[i] = 0;
}

mat_float transpose(mat_float &A) {
    const int n = A.NumRows();
    const int m = A.NumCols();
    mat_float R;
    R.SetDims(m, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            R[j][i] = A[i][j];
    return R;
}

vec_float operator*(double a, const vec_float &b) {
    const int m = b.length();
    vec_float reps;
    reps.SetLength(m);
    clear(reps);
    for (int j = 0; j < m; j++)
        reps[j] = a * b[j];
    return reps;
}


mat_float operator*(const mat_float &A, const mat_float &B) {
    const int n = A.NumRows();
    const int l = A.NumCols();
    const int m = B.NumCols();
    assert(B.NumRows() == l);
    mat_float reps;
    reps.SetDims(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) {
            reps[i][j] = 0;
            for (int k = 0; k < l; k++)
                reps[i][j] += A[i][k] * B[k][j];
        }
    return reps;
}

float operator*(const vec_float &A, const vec_float &B) {
    const int n = A.length();
    assert(B.length() == n);
    float reps = 0;
    for (int i = 0; i < n; i++)
        reps += A[i] * B[i];
    return reps;
}

vec_float operator*(const mat_float &A, const vec_float &b) {
    const int n = A.NumRows();
    const int m = A.NumCols();
    assert(b.length() == m);
    vec_float reps;
    reps.SetLength(n);
    clear(reps);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            reps[i] += A[i][j] * b[j];
    return reps;
}

vec_float operator*(const vec_float &b, const mat_float &A) {
    const int n = A.NumRows();
    const int m = A.NumCols();
    assert(b.length() == n);
    vec_float reps;
    reps.SetLength(m);
    clear(reps);
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
            reps[i] += b[j] * A[j][i];
    return reps;
}

/** diag(w) * A */
mat_float diagProd(const vec_float &w, const mat_float &A) {
    const int n = A.NumRows();
    const int m = A.NumCols();
    assert(w.length() == n);
    mat_float reps;
    reps.SetDims(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            reps[i][j] = w[i] * A[i][j];
    return reps;
}

vec_float operator-(const vec_float &a, const vec_float &b) {
    const int m = b.length();
    assert(a.length() == m);
    vec_float reps;
    reps.SetLength(m);
    clear(reps);
    for (int j = 0; j < m; j++)
        reps[j] = a[j] - b[j];
    return reps;
}

vec_float operator+(const vec_float &a, const vec_float &b) {
    const int m = b.length();
    assert(a.length() == m);
    vec_float reps;
    reps.SetLength(m);
    clear(reps);
    for (int j = 0; j < m; j++)
        reps[j] = a[j] + b[j];
    return reps;
}


/**
 * en entree:
 * G = [ Id | 0 ]
 *     [ 0  | A ]
 * - Id of size i
 * - A is symmetric
 * - alpha = 1/sqrt(G[i][i])
 *
 * this function multiplies row i and column i by alpha.
 * */
void symetricScale(mat_float &G, int i, float alpha) {
    const int k = G.NumRows();
    for (int j = i + 1; j < k; j++) {
        G[i][j] *= alpha;
        G[j][i] = G[i][j];
    }
    G[i][i] = 1.;
}

/**
 * en entree:
 * G = [ Id | 0  ----------- 0 ]
 *     [ 0  | 1 0 -gamma * * * ]
 *     [ 0  | 0 *  * * * * * * ]
 *     ...
 * - Id of size i
 * - G[i][j] = -gamma
 *
 * this function adds gamma times row i to row j.
 * then adds gamma times col i to col j.
 */
void symetricTransvect(mat_float &G, int j, int i, float gamma) {
    const int k = G.NumRows();
    G[j][i] = 0.;
    for (int l = i + 1; l < k; l++) {
        G[j][l] += gamma * G[i][l];
    }
    //just copy the column
    for (int l = i; l < k; l++) {
        G[l][j] = G[j][l];
    }
}


/** multiply row i by alpha */
void rowScale(mat_float &A, int i, float alpha) {
    const int m = A.NumCols();
    for (int j = 0; j < m; j++) A[i][j] *= alpha;
}

/** add gamma times row i to row j */
void rowTransvect(mat_float &A, int j, int i, float gamma) {
    const int m = A.NumCols();
    for (int l = 0; l < m; l++) A[j][l] += gamma * A[i][l];
}

float sqr(float x) { return x * x; }

/** square norm of columns */
vec_float colSqNorms(const mat_float &A) {
    const int n = A.NumRows();
    const int m = A.NumCols();
    vec_float reps;
    reps.SetLength(m);
    clear(reps);
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            reps[j] += sqr(A[i][j]);
        }
    }
    return reps;
}

/** square norm of columns */
vec_float scaledColSqNorms(const mat_float &A, const vec_float &w) {
    const int n = A.NumRows();
    const int m = A.NumCols();
    vec_float reps;
    reps.SetLength(m);
    clear(reps);
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            reps[j] += w[i] * sqr(A[i][j]);
        }
    }
    return reps;
}

float pvalexp_exact(float x) {
    return 1. + erf(-exp(x / 2.) / sqrt(2.));
}

float pvalexp_approx(float x) {
    x = exp(x/2.) / sqrt(2.);
    float a1=.278393;
    float a2=.230389;
    float a3=.000972;
    float a4=.078108;

    float t = (1 + a1 * x + a2 * x*x + a3 * x*x*x + a4 * x*x*x*x);
    t = t*t*t*t;

    return 1.0/t;
}

#define pvalexp pvalexp_exact

void draw_histogram(const string &name, const vec_float &values) {
    static const int NB = 100;
    double vmin = INFINITY;
    double vmax = -INFINITY;
    double vmean = 0;
    double vstdev = 0;
    double n = values.length();
    double bins[NB + 1];
    for (int i = 0; i <= NB; i++) bins[i] = 0;
    for (float v : values) {
        if (v > vmax) vmax = v;
        if (v < vmin) vmin = v;
        vmean += v;
        vstdev += v * v;
    }
    vmean /= n;
    vstdev /= n;
    vstdev -= vmean * vmean;
    vstdev = sqrt(vstdev);
    cout << name << " -- " << vmean << " +- " << vstdev << " -- [ " << vmin << " , " << vmax << " ]" << endl;
    for (float v : values) {
        double x = NB * (v - vmin) / (vmax - vmin);
        double xf = floor(x);
        long xidx = long(xf);
        if (xidx == NB) {
            bins[NB] += 1.;
        } else {
            bins[xidx] += xf + 1 - x;
            bins[xidx + 1] += x - xf;
        }
    }
    ofstream ofs(name + string(".histo"));
    for (int i = 0; i <= NB; i++)
        ofs << vmin + double(i) / double(NB) * (vmax - vmin) << " " << bins[i] << endl;
    ofs.close();
    ofstream ofs2(name + string(".histo.plot"));
    ofs2 << "set term post eps color enhanced" << endl;
    ofs2 << "set output '" << name << ".eps'" << endl;
    ofs2 << "plot '" << name << ".histo' with boxes" << endl;
    ofs2.close();
    system((string("gnuplot ")+name+string(".histo.plot")).c_str());
    system((string("epstopdf ")+name+string(".eps")).c_str());
}

void draw_histogram(const string &name, const mat_float &values) {
    long n = values.NumRows();
    long m = values.NumCols();
    vec_float tmp;
    tmp.SetLength(n * m);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            tmp[i * m + j] = values[i][j];
        }
    }
    draw_histogram(name, tmp);
}


int main() {
    int ITERS = 7;    //num of logreg iters
    double step = 4.; //learning rate (close to 4)

    Data data;
    fill_data(data);

    const int k = data.k;
    const int m = data.m;
    const int n = data.n;
    mat_float& X = data.X;
    mat_float& S = data.S;
    vec_float& y = data.y;

    //
    vec_float p;

    //logistic regression
    vec_float beta;
    beta.SetLength(k);
    clear(beta);
    for (int iter = 0; iter <= ITERS; iter++) {
        cout << "--- iter " << iter << "-----" << endl;
        cout << "beta: " << beta << endl;
        vec_float XBeta = X*beta;
        cout << "XBeta: " << XBeta << endl;
        draw_histogram(string("X_Beta_")+to_string(iter), XBeta);
        p = sigmoid_vec(X * beta);
        cout << "p: " << p << endl;
        vec_float ymp = y-p;
        cout << "ymp: " << ymp << endl;
        vec_float mgrad = (transpose(X) * ymp);
        cout << "mgrad: " << mgrad << endl;
        beta = beta + step * mgrad;
        cout << "grad: " << iter << " " << sqrt(mgrad * mgrad) << endl;
        //cout << "p: " << p << endl;
    }
    // p and W
    vec_float w;
    w.SetLength(n);
    for (int i = 0; i < n; i++) w[i] = p[i] * (1 - p[i]);
    cout << "grad: " << (transpose(X) * (y - p)) << endl;
    cout << "beta: " << beta << endl;
    cout << "p: " << p << endl;
    draw_histogram("p", p);
    cout << "w: " << w << endl;
    draw_histogram("w", w);
    // Zstar
    vec_float zStar;
    zStar = (y - p) * S;
    cout << "zStar:" << zStar << endl;
    draw_histogram("zStar", zStar);
    // G
    mat_float G = transpose(X) * diagProd(w, X);
    mat_float A = transpose(X) * diagProd(w, S);
    cout << "G: " << G << endl;
    draw_histogram("G", G);
    cout << "A: " << A << endl;
    draw_histogram("A", A);
    //cholesky
    for (int i = 0; i < k; i++) {
        assert(sqrt(G[i][i]) > 0); //the matrix must be positive definite!!
        float alpha = 1. / sqrt(G[i][i]);
        symetricScale(G, i, alpha);
        rowScale(A, i, alpha);
        for (int j = i + 1; j < k; j++) {
            float gamma = -G[j][i];
            symetricTransvect(G, j, i, gamma);
            rowTransvect(A, j, i, gamma);
        }
    }
    //denominator
    vec_float sStar2 = scaledColSqNorms(S, w) - colSqNorms(A);
    cout << "sStar2: " << sStar2 << endl;
    draw_histogram("sStar2", sStar2);
    vec_float ri;
    ri.SetLength(m);
    for (int j = 0; j < m; j++)
        ri[j] = 2 * log(abs(zStar[j])) - log(abs(sStar2[j]));
    cout << "ri: " << ri << endl;
    draw_histogram("ri", ri);
    vec_float pval;
    pval.SetLength(m);
    for (int j = 0; j < m; j++)
        pval[j] = pvalexp(ri[j]);
    cout << "pval: " << pval << endl;
    draw_histogram("pval", pval);

    ofstream ofs("pvalexp.dat");
    for (int j = 0; j < m; j++)
        ofs << ri[j] << " " << pval[j] << endl;
    ofs.close();
    return 0;
}
