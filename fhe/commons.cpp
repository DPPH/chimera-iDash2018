#include "commons.h"

NTL_CLIENT;

void ostream_write_binary(std::ostream &out, const void *const data, size_t bytes) {
    out.write((const char *) data, bytes);
}

void istream_read_binary(std::istream &in, void *const data, size_t bytes) {
    in.read((char *) data, bytes);
}

