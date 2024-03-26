#include <bioparser/fasta_parser.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
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

TEST(MyTest, MyFirstTest) { EXPECT_EQ(1, 1); }

TEST(MyTest, MySecondTest) {
  const char *triobinningPath = getenv("TRIOBINNING_PATH");
  if (!triobinningPath) {
    FAIL() << "TRIOBINNING_PATH environment variable is not set";
    return;
  }

  string pathRef = string(triobinningPath) + "ecoli_reads.fastq";
  auto p1 =
      bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(pathRef);
  auto ref = p1->Parse(-1);
  EXPECT_EQ(ref.size(), 2);
}