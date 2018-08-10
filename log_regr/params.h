#ifndef PARAMS_H
#define PARAMS_H

#include <cmath>
#include <stdio.h>

struct LRParams
{
    // data set shape
    int n;          // number X rows
    int k;          // number X cols
    int m;          // number SNPs

    int X_precbit = 8;

    //
    int nb_iters;
    int nb_threads;
    int seed;
    int verbose_level;

    int sigm_levels = 16;
    double sigm_in_min = -2.0;
    double sigm_in_max = 2.0;
    double scale = 1.0 / 2;

    const char* const params_filename = "params.bin";
    const char* const secret_keyset_filename = "secret_keyset.bin";
    const char* const cloud_keyset_filename = "cloud_keyset.bin";

    /**
     * @brief Update computed options
     */
    void update() {
    }

    void print() {
        printf("LRParams\n");

        printf("  data:\n");
        printf("    %d individuals\n", n);
        printf("    %d features\n", k-1);
        printf("    %d SNPs\n", m);

        printf("  algorithm:\n");
        printf("    %d iterations\n", nb_iters);
        printf("    %d threads\n", nb_threads);
        printf("    %d seed\n", seed);
        printf("    %d verbose level\n", verbose_level);
    }
};

#endif
