#include "kmer_maker.hpp"

#include <biosoup/nucleic_acid.hpp>
#include <cstdint>
#include <future>
#include <iostream>
#include <iterator>
#include <math.h>
#include <string>
#include <thread_pool/thread_pool.hpp>
#include <unordered_map>
#include <vector>
#include <set>

using namespace std;

namespace biosoup {
atomic<uint32_t> NucleicAcid::num_objects{0};
} // namespace biosoup

namespace kmer {

uint64_t invertibleHash(uint64_t key) {

  // Pseudo-code from https://naml.us/post/inverse-of-a-hash-function/
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}

unsigned int returnHash(vector<uint64_t> list, unsigned int kmer_length) {
  // Concatenate the numbers in the list
  string concatenated_number = "";
  for (const auto &num : list) {
    concatenated_number += to_string(num);
  }
  // Return the hash of the concatenated number
  return invertibleHash(stoull(concatenated_number));
}

set<unsigned int> kmer_maker_set(std::string genome, unsigned int kmer_length) {
  // Make a set of kmers
  set<unsigned int> kmer_set;

  // Iterate through the genome
  for (unsigned int i = 0; i < genome.size() - kmer_length + 1; i++) {
    // Get the kmer
    string kmer = genome.substr(i, kmer_length);

    // Get the hash of the kmer
    biosoup::NucleicAcid na("name", kmer);
    auto list = na.deflated_data;
    unsigned int hash = returnHash(list, kmer_length);

    // Add the hash to the set
    kmer_set.insert(hash);
  }
  return kmer_set;
}

unordered_map<unsigned int, unsigned int>
kmer_maker(vector<unique_ptr<seq::Sequence>> &ref, unsigned int kmer_length) {

  // Make a map of kmers and their counts
  unordered_map<unsigned int, unsigned int> kmer_map;

  // Iterate through the reference sequences
  for (auto &seq : ref) {
    // Iterate through the sequence
    for (unsigned int i = 0; i < seq->genome.size() - kmer_length + 1; i++) {

      // Get the kmer
      string kmer = seq->genome.substr(i, kmer_length);

      // Get the hash of the kmer
      biosoup::NucleicAcid na(seq->name, kmer);
      auto list = na.deflated_data;
      unsigned int hash = returnHash(list, kmer_length);

      // Add the kmer to the map
      if (kmer_map.find(hash) == kmer_map.end()) {
        kmer_map[hash] = 1;
      } else {
        kmer_map[hash]++;
      }
    }
  }
  return kmer_map;
}

unordered_map<unsigned int, unsigned int>
parallel_kmer(vector<unique_ptr<seq::Sequence>> &ref, unsigned int kmer_length,
              unsigned int num_threads) {

  // Initialize thread_pool
  vector<future<unordered_map<unsigned int, unsigned int>>> thread_features;
  shared_ptr<thread_pool::ThreadPool> thread_pool_ =
      make_shared<thread_pool::ThreadPool>(num_threads);

  // Split the reference sequences into parts
  const size_t part_size = ref.size() / num_threads;
  for (uint32_t i = 0; i < num_threads; ++i) {
    thread_features.emplace_back(thread_pool_->Submit(
        [&ref, kmer_length, i, num_threads, part_size]() {

          // Get the part of the reference sequences
          size_t start_idx = i * part_size;
          size_t end_idx = (i == num_threads - 1) ? ref.size() : (start_idx + part_size);
          vector<unique_ptr<seq::Sequence>> part_ref(
              make_move_iterator(ref.begin() + start_idx),
              make_move_iterator(ref.begin() + end_idx));

          // Use kmer_maker to get the kmers for the part of the reference sequences
          auto kmer_maps = kmer_maker(part_ref, kmer_length);
          return kmer_maps;
        }));
  }

  // Initialize the map of kmers and their counts
  unordered_map<unsigned int, unsigned int> kmer_map = {};

  // Combine the results from the threads
  for (auto &thread_feature : thread_features) {
    auto map = thread_feature.get();
    if (kmer_map.empty()) {
      kmer_map = map;
      continue;
    }
    for (auto &pair : map) {
      if (kmer_map.find(pair.first) == kmer_map.end()) {
        kmer_map[pair.first] = pair.second;
      } else {
        kmer_map[pair.first] += pair.second;
      }
    }
  }
  thread_features.clear();

  // Return the map of kmers and their counts
  return kmer_map;
}
} // namespace kmer