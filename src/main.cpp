#include "sequence.hpp"
#include <algorithm>
#include <bioparser/fasta_parser.hpp>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <kmer_finder.hpp>
#include <kmer_maker.hpp>
#include <memory>
#include <set>
#include <sstream>
#include <unordered_map>
#include <vector>

using namespace std;
const string path =
    getenv("TRIOBINNING_PATH") ? getenv("TRIOBINNING_PATH") : ".";
const int KMER_LENGTH = 16;
const int NUM_THREADS = 6;

int main() {
  string pathRef = path + "external/genomes/ecoli_reads.fastq";
  string pathRef2 = path + "external/genomes/ecoli_mutated_reads.fastq";
  string pathRef3 = path + "external/genomes/output_reads.fastq";
  auto p1 =
      bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(pathRef);
  auto ref = p1->Parse(-1);

  auto p2 = bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(
      pathRef2);
  auto ref2 = p2->Parse(-1);

  auto p3 = bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(
      pathRef3);
  auto ref3 = p3->Parse(-1);

  auto start = chrono::high_resolution_clock::now();

  auto map = kmer::parallel_kmer(ref, KMER_LENGTH, NUM_THREADS);
  auto map2 = kmer::parallel_kmer(ref2, KMER_LENGTH, NUM_THREADS);

  auto end = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::seconds>(end - start).count();

  cout << "Time taken: " << duration << " seconds" << endl;
  cout << "Number of kmers: " << map.size() << endl;
  cout << "Number of kmers: " << map2.size() << endl;

  auto currentTime = std::chrono::high_resolution_clock::now();
  auto currentTime_t = std::chrono::system_clock::to_time_t(currentTime);
  std::cout << "Start time: " << std::ctime(&currentTime_t) << std::endl;

  set<int> keys1;
  set<int> keys2;

  for (auto &it : map) {
    keys1.insert(it.first);
  }

  for (auto &it : map2) {
    keys2.insert(it.first);
  }

  set<int> diff;
  set_difference(keys1.begin(), keys1.end(), keys2.begin(), keys2.end(),
                 inserter(diff, diff.begin()));

  set<int> diff2;
  set_difference(keys2.begin(), keys2.end(), keys1.begin(), keys1.end(),
                 inserter(diff2, diff2.begin()));

  map.clear();
  map2.clear();

  auto myMap =
      kmer::parallel_finder(ref3, KMER_LENGTH, NUM_THREADS, keys1, keys2);

  cout << diff.size() << endl << diff2.size() << endl;

  int count = 0;
  for (const auto &pair : myMap) {
    std::cout << pair.first << ": " << pair.second << std::endl;
    count++;
    if (count == 10) {
      break; // Stop printing after 10 elements
    }
  }
  /*
  for (auto &item : map) {
    if (diff.find(item.first) == diff.end()) {
      if (map[item.first] > map2[item.first])
        diff.insert(item.first);
      else
        diff2.insert(item.first);
    }
  }
  */

  auto end2 = chrono::high_resolution_clock::now();
  auto duration2 = chrono::duration_cast<chrono::seconds>(end2 - end).count();

  cout << "Time taken to make differences: " << duration2 << " seconds" << endl;
  // cout << "Number of differences: " << diff.size() << endl;
  // cout << "Number of differences: " << diff2.size() << endl;

  return 0;
}
