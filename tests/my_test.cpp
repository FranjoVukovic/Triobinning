#include "kmer_maker.hpp"
#include "kmer_finder.hpp"
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

  string pathRef = string(triobinningPath) + "ecoli_reads.fasta";
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

  string pathRef = string(triobinningPath) + "ecoli_reads.fasta";
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

  string pathRef = string(triobinningPath) + "ecoli_reads.fasta";
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

  string pathRef = string(triobinningPath) + "ecoli_reads.fasta";
  auto p1 =
      bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(pathRef);
  auto ref = p1->Parse(-1);

  auto map = kmer::parallel_kmer(ref, KMER_LENGTH, 4);

  EXPECT_EQ(map.size(), 32200);
}

TEST(Classification, MultiThreded) {
  const char *triobinningPath = getenv("TRIOBINNING_PATH");
  if (!triobinningPath) {
    FAIL() << "TRIOBINNING_PATH environment variable is not set";
    return;
  }

  string pathRef = string(triobinningPath) + "ecoli_reads.fasta";
  string pathRef2 = string(triobinningPath) + "ecoli_mutated_reads.fasta";
  string pathRef3 = string(triobinningPath) + "output.fasta";

  auto p1 =
      bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(pathRef);
  auto ref = p1->Parse(-1);

  auto map = kmer::parallel_kmer(ref, KMER_LENGTH, 4);
  set<int> keys1;
  for (auto &pair : map) {
    keys1.insert(pair.first);
  }
  map.clear();
  ref.clear();

  auto p2 =
      bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(pathRef2);
  auto ref2 = p2->Parse(-1);

  auto map2 = kmer::parallel_kmer(ref2, KMER_LENGTH, 4);
  set<int> keys2;
  for (auto &pair : map2) {
    keys2.insert(pair.first);
  }
  map2.clear();
  ref2.clear();

  set<int> diff;
  set_difference(keys1.begin(), keys1.end(), keys2.begin(), keys2.end(),
                 inserter(diff, diff.begin()));

  set<int> diff2;
  set_difference(keys2.begin(), keys2.end(), keys1.begin(), keys1.end(),
                 inserter(diff2, diff2.begin()));

  keys1.clear();
  keys2.clear();

  auto p3 =
      bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(pathRef3);
  auto ref3 = p3->Parse(-1);

  auto myMap = kmer::parallel_finder(ref3, KMER_LENGTH, 4, diff, diff2);

  int paternal = 0;
  int maternal = 0;
  int unknown = 0;

  for (auto &pair : myMap) {
    if (pair.second == 1) {
      paternal++;
    } else if (pair.second == 2) {
      maternal++;
    } else {
      unknown++;
    }
  }

  EXPECT_EQ(paternal, 3);
  EXPECT_EQ(maternal, 3);
  EXPECT_EQ(unknown, 0);
}