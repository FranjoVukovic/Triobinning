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

using namespace std;

namespace biosoup {
atomic<uint32_t> NucleicAcid::num_objects{0};
} // namespace biosoup

namespace kmer {

uint64_t invertibleHash(uint64_t key) {
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
  string concatenated_number = "";
  for (const auto &num : list) {
    concatenated_number += to_string(num);
  }
  return invertibleHash(stoull(concatenated_number));
}

unordered_map<unsigned int, unsigned int>
kmer_maker(vector<unique_ptr<seq::Sequence>> &ref, unsigned int kmer_length) {

  unordered_map<unsigned int, unsigned int> kmer_map;
  for (auto &seq : ref) {
    for (unsigned int i = 0; i < seq->genome.size() - kmer_length + 1; i++) {
      string kmer = seq->genome.substr(i, kmer_length);
      biosoup::NucleicAcid na(seq->name, kmer);
      auto list = na.deflated_data;
      unsigned int hash = returnHash(list, kmer_length);
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

  vector<future<unordered_map<unsigned int, unsigned int>>> thread_features;
  shared_ptr<thread_pool::ThreadPool> thread_pool_ =
      make_shared<thread_pool::ThreadPool>(num_threads);

  const size_t part_size = ref.size() / num_threads;

  vector<unique_ptr<seq::Sequence>> ref_copy;
  ref_copy.reserve(ref.size());
  for (const auto &seq : ref) {
    ref_copy.push_back(make_unique<seq::Sequence>(*seq));
  }

  for (uint32_t i = 0; i < num_threads; ++i) {
    thread_features.emplace_back(thread_pool_->Submit(
        [&ref_copy, kmer_length, i, num_threads, part_size]() {
          size_t start_idx = i * part_size;
          size_t end_idx = (i == num_threads - 1) ? ref_copy.size()
                                                  : (start_idx + part_size);
          unordered_map<unsigned int, unsigned int> kmer_maps;
          vector<unique_ptr<seq::Sequence>> part_ref(
              make_move_iterator(ref_copy.begin() + start_idx),
              make_move_iterator(ref_copy.begin() + end_idx));

          for (auto &seq : part_ref) {
            for (unsigned int i = 0; i < seq->genome.size() - kmer_length + 1;
                 i++) {
              string kmer = seq->genome.substr(i, kmer_length);
              biosoup::NucleicAcid na(seq->name, kmer);
              auto list = na.deflated_data;
              unsigned int hash = returnHash(list, kmer_length);
              if (kmer_maps.find(hash) == kmer_maps.end()) {
                kmer_maps[hash] = 1;
              } else {
                kmer_maps[hash]++;
              }
            }
          }
          return kmer_maps;
        }));
  }

  unordered_map<unsigned int, unsigned int> kmer_map = {};

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

  return kmer_map;
}
} // namespace kmer