#include "commons.h"
#include <map>

NTL_CLIENT;

void ostream_write_binary(std::ostream &out, const void *const data, size_t bytes) {
    out.write((const char *) data, bytes);
}

void istream_read_binary(std::istream &in, void *const data, size_t bytes) {
    in.read((char *) data, bytes);
}

void store_forever(std::shared_ptr<void> object) {
    static std::map<UINT64, std::shared_ptr<void>> v;
    UINT64 addr = (UINT64) object.get();
    if (v.find(addr) == v.end()) {
        v.emplace(addr, object);
    }
}

