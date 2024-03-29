#include "kmer_maker.hpp"

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

namespace kmer {

const int BIT_LENGTH = 32;

unsigned int invertibleHash(char c, int p) {
  uint32_t x = static_cast<uint32_t>(c); // Convert char to uint32_t
  uint32_t m = (std::pow(2, p)) - 1;
  x = (~x + (x << 21)) & m;
  x = x ^ (x >> 24);
  x = (x + (x << 3) + (x << 8)) & m;
  x = x ^ (x >> 14);
  x = (x + (x << 2) + (x << 4)) & m;
  x = x ^ (x >> 28);
  x = (x + (x << 31)) & m;
  return x;
}

unsigned int returnHash(const string &kmer, unsigned int kmer_length) {
  unsigned int hash = 0;
  for (unsigned int i = 0; i < kmer_length; i++) {
    hash += invertibleHash(kmer[i], BIT_LENGTH);
  }
  return hash;
}

unordered_map<unsigned int, unsigned int>
kmer_maker(vector<unique_ptr<seq::Sequence>> &ref, unsigned int kmer_length) {

  unordered_map<unsigned int, unsigned int> kmer_map;
  for (auto &seq : ref) {
    for (unsigned int i = 0; i < seq->genome.size() - kmer_length + 1; i++) {
      string kmer = seq->genome.substr(i, kmer_length);
      unsigned int hash = returnHash(kmer, kmer_length);
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

  unordered_map<unsigned int, unsigned int> kmer_map = {};
  vector<future<unordered_map<unsigned int, unsigned int>>> thread_features;
  shared_ptr<thread_pool::ThreadPool> thread_pool_ =
      make_shared<thread_pool::ThreadPool>(num_threads);

  const size_t part_size = ref.size() / num_threads;

  for (uint32_t i = 0; i < num_threads; ++i) {
    thread_features.emplace_back(
        thread_pool_->Submit([&ref, kmer_length, i, num_threads, part_size]() {
          size_t start_idx = i * part_size;
          size_t end_idx =
              (i == num_threads - 1) ? ref.size() : (start_idx + part_size);
          unordered_map<unsigned int, unsigned int> kmer_maps;
          vector<unique_ptr<seq::Sequence>> part_ref(
              make_move_iterator(ref.begin() + start_idx),
              make_move_iterator(ref.begin() + end_idx));
          for (auto &seq : part_ref) {
            for (unsigned int i = 0; i < seq->genome.size() - kmer_length + 1;
                 i++) {
              string kmer = seq->genome.substr(i, kmer_length);
              unsigned int hash = returnHash(kmer, kmer_length);
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