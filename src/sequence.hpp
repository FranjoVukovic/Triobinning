#pragma once

#include <cstdint>
#include <string>

using namespace std;

namespace seq {
struct Sequence {
public:
  Sequence(const char *name, uint32_t name_len, const char *genome,
           uint32_t gen_len) {
    this->name = string(name, name_len);
    this->genome = string(genome, gen_len);
  }

  string name;
  string genome;
};
} // namespace seq
