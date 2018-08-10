#ifndef FHE_COMMONS_H
#define FHE_COMMONS_H

#include <NTL/LLL.h>

#if __APPLE__
typedef unsigned long long UINT64;
#else
typedef uint64_t UINT64;
#endif

#define NA 0x8000000000000000ul
const int64_t BITS_PER_LIMBS = 64;
//const int64_t BYTES_PER_LIMBS = BITS_PER_LIMBS / 8;

inline int64_t limb_precision(UINT64 bit_precision) {
    if (bit_precision == NA) return NA;
    return (bit_precision + BITS_PER_LIMBS - 1) / BITS_PER_LIMBS;
}

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

#define NO_COPY(TypeName) \
    TypeName(const TypeName&)=delete; \
    void operator=(const TypeName&)=delete

#endif //FHE_COMMONS_H
