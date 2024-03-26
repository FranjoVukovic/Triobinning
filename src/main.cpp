#include <algorithm>
#include <bioparser/fasta_parser.hpp>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

using namespace std;
const string path =
    getenv("TRIOBINNING_PATH") ? getenv("TRIOBINNING_PATH") : ".";

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

int main() {
  string pathRef = path + "ecoli.fna";
  auto p1 =
      bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(pathRef);
  auto ref = p1->Parse(-1);
  const Sequence &refrence = *ref[0];

  cout << refrence.name << endl;

  return 0;
}