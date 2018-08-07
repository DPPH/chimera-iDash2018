#ifndef IO_H
#define IO_H

#include <iostream>
#include <NTL/LLL.h>
#include <cassert>
#include <fstream>
#include <sstream>

NTL_CLIENT;

typedef NTL::Vec<float> vec_float;
typedef NTL::Mat<float> mat_float;

using namespace std;

struct Data {
  int n;        // number X rows
  int k;        // number X cols
  int m;        // number SNPs

  mat_float X;
  vec_float y;
  mat_float S;
};

bool is_binary(float x) {
    return x == 0 || x == 1;
}

void read_S(string filename, Data& data) {
    ifstream ifs(filename); //"data/snpMat.txt"
    assert(ifs);

    mat_float &S = data.S;
    const int n = data.n;
    const int m = data.m;

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

void read_ortho_Xy(std::string filename, Data& data) {
    ifstream ifs(filename); // "data/covariates.csv"
    assert(ifs);

    mat_float &X = data.X;
    vec_float &y = data.y;
    const int n = data.n;
    const int k = data.k;

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

void fill_data(Data& data) {
  data.n = 245;
  data.k = 4;
  data.m = 10643;

  read_S("data/snpMat.txt", data);
  read_ortho_Xy("data/covariates.csv", data);
}

#endif
