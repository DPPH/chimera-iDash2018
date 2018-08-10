#ifndef COMMON_H
#define COMMON_H

#include <cmath>

float sigmoid(float x) {
    return 1. / (1 + exp(-x));
}

#endif
