#include <iostream>
#include <NTL/LLL.h>
#include <NTL/ZZ_limbs.h>

NTL_CLIENT;

typedef NTL::Vec<float> vec_float;
typedef NTL::Mat<float> mat_float;


float sigmoid(float x) {
    return 1./(1+exp(-x));
}

vec_float sigmoid_vec(const vec_float& x) {
    const int n = x.length();
    vec_float reps; reps.SetLength(n);
    for (int i=0; i<n; i++) reps[i]=sigmoid(x[i]);
}

int main() {
    int k = 4;
    int m = 10000;
    int n = 235;
    int ITERS = 5;    //num of logreg iters
    double step = 5.; //learning rate
    mat_float X;
    vec_float y;
    //
    vec_float p;

    NTL::Vec<float> beta; beta.SetLength(n); clear(beta);
    for (int iter=0; iter<=ITERS; iter++) {
        p =
        beta = beta + step * transpose(X) * (y - sigmoid_vec(X * beta));
    }
    vec_float p =
    mat_ZZ m;
    m.SetDims(3,4);
    clear(m);
    m(1,1)=123;
    cout << ZZ_limbs_get(m(1,1))[0] << endl;
    std::cout << "Hello, World!"  << m << std::endl;
    return 0;
}