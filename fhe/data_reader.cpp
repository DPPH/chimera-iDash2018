#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include "section2_params.h"

using namespace std;

int64_t read_dims_from_data(int64_t index) {
    static int64_t dims[4]; // n, k, m, lvl0_n
    static bool initialized = false;
    if (!initialized) {
        //try to initialize the dimensions from the dataset if it is available, else set it from the dimensions.txt
        {
            ifstream ifs("data/covariates.csv");
            if (!ifs) {
                //try to read from dimensions.txt
                ifstream dfs("dimensions.txt");
                assert_dramatically(dfs, "dimensions.txt not found");
                string buf;
                dfs >> buf >> dims[0]
                    >> buf >> dims[1]
                    >> buf >> dims[2]
                    >> buf >> dims[3];
                dfs.close();
                cout << "dimensions: n,k,m,n_lvl0: " << dims[0] << "," << dims[1]
                     << "," << dims[2] << "," << dims[3] << endl;
                initialized = true;
                return dims[index];
            }
        }
        //read n and k from covariates.csv
        string line;
        int64_t n = -1;
        int64_t k = -1;
        int64_t m = 0;
        int64_t n_lvl0 = 0;
        {
            ifstream ifs("data/covariates.csv");
            assert_dramatically(ifs, "data/covariates.csv not found");
            for (std::getline(ifs, line); ifs; std::getline(ifs, line)) {
                n++; //count number of lines
                if (n == 1) {
                    //count number of ,
                    for (char c: line) {
                        if (c == ',') k++;
                    }
                }
            }
            ifs.close();
        }
        //read m
        {
            ifstream ifs2("data/snpMat.txt");
            assert_dramatically(ifs2, "data/covariates.csv not found");
            std::getline(ifs2, line);
            std::getline(ifs2, line);
            int64_t buf = 0;
            istringstream iss(line);
            for (iss >> buf; iss; iss >> buf) {
                assert_dramatically(buf == 0 || buf == 1, "snpMat is not binary");
                m++;
            }
            ifs2.close();
        }
        //read lvl0_n
        {
            ifstream ifs3("params.bin");
            assert_dramatically(ifs3, "params.bin not found (have you run section 1 keygen?)");
            std::getline(ifs3, line);
            std::getline(ifs3, line);
            std::getline(ifs3, line);
            std::getline(ifs3, line);
            string tmp;
            istringstream iss(line);
            iss >> tmp >> n_lvl0;
            assert_dramatically(tmp == string("n:"));
        }
        dims[0] = n;
        dims[1] = k;
        dims[2] = m;
        dims[3] = n_lvl0;
        ofstream ofs("dimensions.txt");
        ofs << "n: " << n << endl
            << "k: " << k << endl
            << "m: " << m << endl
            << "n_lvl0: " << n_lvl0 << endl;
        ofs.close();
        cout << "dimensions: n,k,m,n_lvl0: " << n << "," << k << "," << m << "," << n_lvl0 << endl;
        initialized = true;
    }
    return dims[index];
}


