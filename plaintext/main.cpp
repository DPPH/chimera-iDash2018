#include <iostream>
#include <NTL/LLL.h>
#include <NTL/ZZ_limbs.h>
#include <cassert>
#include <fstream>
#include <sstream>

NTL_CLIENT;

typedef NTL::Vec<float> vec_float;
typedef NTL::Mat<float> mat_float;


float sigmoid(float x) {
    return 1. / (1 + exp(-x));
}

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
    return 1. + erf(-exp(x/2.) / sqrt(2.));
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

bool is_binary(float x) {
    return x == 0 || x == 1;
}

void fill_matrix_S(mat_float &S, int n, int m) {
    ifstream ifs("data/snpMat.txt");
    assert(ifs);
    S.SetDims(n, m);
    string line;
    std::getline(ifs, line); //ignore first header line
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) {
            ifs >> S[i][j];
            assert(is_binary(S[i][j])); //S binary
        }
    assert(ifs); //verify that all values have been read
    ifs.close();
}

void fill_matrix_Xy(mat_float &X, vec_float &y, int n, int k) {
    ifstream ifs("data/covariates.csv");
    assert(ifs);
    X.SetDims(n, k);
    y.SetLength(n);
    mat_RR B;
    vec_RR sums;
    vec_RR nbels;
    B.SetDims(k, n);
    sums.SetLength(k);
    clear(sums);
    nbels.SetLength(k);
    clear(nbels);

    string line;
    string buf;
    std::getline(ifs, line); //ignore first header line
    for (int i = 0; i < n; i++) {
        std::getline(ifs, line);
        for (int j = 0; j < int(line.size()); j++) {
            if (line[j] == ',') line[j] = ' ';
        }
        istringstream iss(line);
        iss >> buf; //ignore label
        iss >> y[i]; //read y
        assert(is_binary(y[i])); //S binary
        B[0][i] = 1.; //intercept first
        nbels[0]++;
        for (int j = 1; j < k; j++) {
            iss >> buf;
            if (buf == "NA") {
                B[j][i] = -1e80;
            } else {
                B[j][i] = stod(buf);
                sums[j] += B[j][i];
                nbels[j]++;
            }
        }
        assert(iss);
    }
    assert(ifs);
    ifs.close();
    //replace all NANs by average
    //cout << B << endl;
    //cout << nbels << endl;
    for (int j = 0; j < k; j++) {
        sums[j] /= nbels[j];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            if (B[j][i] < 1e-75) {
                B[j][i] = sums[j];
            }
        }
    }
    //orthogonalize B
    for (int i = 0; i < k; i++) {
        //remove component on previous vectors
        for (int j = 0; j < i; j++) {
            B[i] -= (B[i] * B[j]) * B[j];
        }
        //normalize B[i]
        RR alpha = inv(sqrt(B[i] * B[i]));
        B[i] *= alpha;
    }
    //cout << B << endl;
    //cout << B*transpose(B) << endl;
    //put B in X
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            conv(X[i][j], B[j][i]);
        }
    }
}


int main() {
    int k = 4;
    int m = 10643;
    int n = 245;
    int ITERS = 10;    //num of logreg iters
    double step = 5.; //learning rate
    mat_float X;
    mat_float S;
    vec_float y;
    fill_matrix_S(S, n, m);
    fill_matrix_Xy(X, y, n, k);

    //
    vec_float p;

    //logistic regression
    vec_float beta;
    beta.SetLength(k);
    clear(beta);
    for (int iter = 0; iter <= ITERS; iter++) {
        p = sigmoid_vec(X * beta);
        vec_float mgrad = (transpose(X) * (y - p));
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
    cout << "w: " << w << endl;
    // Zstar
    vec_float zStar;
    zStar = (y - p) * S;
    cout << "zStar:" << zStar << endl;
    // G
    mat_float G = transpose(X) * diagProd(w, X);
    mat_float A = transpose(X) * diagProd(w, S);
    cout << "G: " << G << endl;
    cout << "A: " << A << endl;
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
    vec_float ri;
    ri.SetLength(m);
    for (int j = 0; j < m; j++)
        ri[j] = 2 * log(abs(zStar[j])) - log(abs(sStar2[j]));
    cout << "ri: " << ri << endl;
    vec_float pval;
    pval.SetLength(m);
    for (int j = 0; j < m; j++)
        pval[j] = pvalexp(ri[j]);
    cout << "pval: " << pval << endl;

    ofstream ofs("pvalexp.dat");
    for (int j = 0; j < m; j++)
        ofs << ri[j] << " " << pval[j] << endl;
    ofs.close();
    return 0;
}
