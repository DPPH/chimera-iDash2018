#ifndef FHE_COMMONS_H
#define FHE_COMMONS_H

#include <NTL/LLL.h>
#include <memory>
#include <vector>

#if __APPLE__
typedef unsigned long UINT64;
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

/**
 * wrapper to ostream.read, that was wrongly defined with char* instead of void*
 * @param out the output stream to write to
 * @param data the data that needs to be written
 * @param bytes the number of bytes to write
 */
void ostream_write_binary(std::ostream &out, const void *const data, size_t bytes);

/**
 * wrapper to istream.write, that was wrongly defined with char* instead of void*
 * @param in the input stream to read from
 * @param data the location where to put the result
 * @param bytes the number of bytes to read
 */
void istream_read_binary(std::istream &in, void *const data, size_t bytes);

template<typename T>
void store_forever(std::shared_ptr<T> object) {
    static std::vector<std::shared_ptr<void>> v;
    v.push_back(object);
}


#define NO_COPY(TypeName) \
    TypeName(const TypeName&)=delete; \
    void operator=(const TypeName&)=delete


static const int64_t BIGTORUS_PARAMS_SERIAL_ID = 1234567;
static const int64_t BIGTORUS_SERIAL_ID = 1234568;
static const int64_t BIGTORUSVECTOR_SERIAL_ID = 1234569;
static const int64_t TLWE_PARAMS_SERIAL_ID = 1234570;
static const int64_t TLWE_KEY_SERIAL_ID = 1234571;
static const int64_t TLWE_SERIAL_ID = 1234572;
static const int64_t TRLWE_PARAMS_SERIAL_ID = 1234573;
static const int64_t TRLWE_SERIAL_ID = 1234574;
static const int64_t PUBKSKEY32_SERIAL_ID = 1234575;
static const int64_t BIG_REAL_ID = 1234576;
static const int64_t BIG_COMPLEX_ID = 1234577;

#endif //FHE_COMMONS_H
