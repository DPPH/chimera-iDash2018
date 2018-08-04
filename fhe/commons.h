#ifndef FHE_COMMONS_H
#define FHE_COMMONS_H

#include <NTL/LLL.h>

const int64_t NA = 1ul << 63ul;
const int64_t BITS_PER_LIMBS = 64;
const int64_t BYTES_PER_LIMBS = BITS_PER_LIMBS / 8;

inline void assert_dramatically(bool condition, const std::string &message = "") {
    if (!condition) {
        std::cerr << "ASSERT_FAILED: " << message << std::endl;
        abort();
    }
}

inline void assert_weakly(bool condition, const std::string &message = "") {
    if (!condition) {
        std::cerr << "ASSERT_WARNING: " << message << std::endl;
    }
}

#endif //FHE_COMMONS_H
