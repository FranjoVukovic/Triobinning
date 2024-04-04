#include "kmer_maker.hpp"
#include <bioparser/fasta_parser.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

const uint32_t KMER_LENGTH = 16;

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
  EXPECT_EQ(ref.size(), 3);
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

  auto map = kmer::kmer_maker(ref, KMER_LENGTH);

  EXPECT_EQ(map.size(), 32200);
}

TEST(MakingKmers, OnceThreaded) {
  const char *triobinningPath = getenv("TRIOBINNING_PATH");
  if (!triobinningPath) {
    FAIL() << "TRIOBINNING_PATH environment variable is not set";
    return;
  }

  string pathRef = string(triobinningPath) + "ecoli_reads.fastq";
  auto p1 =
      bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(pathRef);
  auto ref = p1->Parse(-1);

  auto map = kmer::parallel_kmer(ref, KMER_LENGTH, 1);

  EXPECT_EQ(map.size(), 32200);
}

TEST(MakingKmers, MultiThreaded) {
  const char *triobinningPath = getenv("TRIOBINNING_PATH");
  if (!triobinningPath) {
    FAIL() << "TRIOBINNING_PATH environment variable is not set";
    return;
  }

  string pathRef = string(triobinningPath) + "ecoli_reads.fastq";
  auto p1 =
      bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(pathRef);
  auto ref = p1->Parse(-1);

  auto map = kmer::parallel_kmer(ref, KMER_LENGTH, 4);

  EXPECT_EQ(map.size(), 32200);
}