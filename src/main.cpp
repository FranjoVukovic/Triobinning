#include "sequence.hpp"
#include <algorithm>
#include <bioparser/fasta_parser.hpp>
#include <chrono>
#include <cstring>
#include <iostream>
#include <kmer_maker.hpp>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <vector>

using namespace std;
const string path =
    getenv("TRIOBINNING_PATH") ? getenv("TRIOBINNING_PATH") : ".";

int main() {
  string pathRef = path + "external/genomes/ecoli_reads.fastq";
  auto p1 =
      bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(pathRef);
  auto ref = p1->Parse(-1);

  auto start = chrono::high_resolution_clock::now();

  auto map = kmer::parallel_kmer(ref, 31, 6);

  auto end = chrono::high_resolution_clock::now();

  auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);

  cout << "Time taken: " << duration.count() / 1000 << " seconds" << endl;

  cout << "Number of kmers: " << map.size() << endl;

  return 0;
}
