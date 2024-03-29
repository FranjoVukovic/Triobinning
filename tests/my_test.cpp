#include "kmer_maker.hpp"
#include <bioparser/fasta_parser.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

TEST(EntryTest, ChekingInput) {
  const char *triobinningPath = getenv("TRIOBINNING_PATH");
  if (!triobinningPath) {
    FAIL() << "TRIOBINNING_PATH environment variable is not set";
    return;
  }

  string pathRef = string(triobinningPath) + "ecoli_reads.fastq";
  auto p1 =
      bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(pathRef);
  auto ref = p1->Parse(-1);
  EXPECT_EQ(ref.size(), 2);
}

TEST(MakingKmers, NonThreaded) {
  const char *triobinningPath = getenv("TRIOBINNING_PATH");
  if (!triobinningPath) {
    FAIL() << "TRIOBINNING_PATH environment variable is not set";
    return;
  }

  string pathRef = string(triobinningPath) + "ecoli_reads.fastq";
  auto p1 =
      bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(pathRef);
  auto ref = p1->Parse(-1);

  auto map = kmer::kmer_maker(ref, 31);

  EXPECT_EQ(map.size(), 20164);
}

TEST(MakingKmers, Threaded) {
  const char *triobinningPath = getenv("TRIOBINNING_PATH");
  if (!triobinningPath) {
    FAIL() << "TRIOBINNING_PATH environment variable is not set";
    return;
  }

  string pathRef = string(triobinningPath) + "ecoli_reads.fastq";
  auto p1 =
      bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(pathRef);
  auto ref = p1->Parse(-1);

  auto map = kmer::parallel_kmer(ref, 31, 1);

  EXPECT_EQ(map.size(), 20164);
}