#include "kmer_finder.hpp"
#include "kmer_maker.hpp"

#include <future>
#include <iostream>
#include <set>
#include <string>
#include <thread_pool/thread_pool.hpp>
#include <unordered_map>
#include <vector>

namespace kmer {

unordered_map<unsigned int, std::string>
parallel_finder(vector<unique_ptr<seq::Sequence>> &ref,
                unsigned int kmer_length, unsigned int num_threads,
                set<int> &difference1, set<int> &difference2) {

  unordered_map<unsigned int, std::string> kmer_map;
  vector<future<unordered_map<unsigned int, std::string>>> thread_features;
  shared_ptr<thread_pool::ThreadPool> thread_pool_ =
      make_shared<thread_pool::ThreadPool>(num_threads);

  const size_t part_size = ref.size() / num_threads;

  vector<unique_ptr<seq::Sequence>> ref_copy;
  ref_copy.reserve(ref.size());
  for (const auto &seq : ref) {
    ref_copy.push_back(make_unique<seq::Sequence>(*seq));
  }

  for (uint32_t i = 0; i < num_threads; ++i) {
    thread_features.emplace_back(
        thread_pool_->Submit([&ref_copy, kmer_length, i, num_threads, part_size,
                              difference1, difference2]() {
          size_t start_idx = i * part_size;
          size_t end_idx = (i == num_threads - 1) ? ref_copy.size()
                                                  : (start_idx + part_size);
          unordered_map<unsigned int, std::string> kmer_maps;
          vector<unique_ptr<seq::Sequence>> part_ref(
              make_move_iterator(ref_copy.begin() + start_idx),
              make_move_iterator(ref_copy.begin() + end_idx));
          auto map = kmer::kmer_maker(part_ref, kmer_length);
          for (auto &it : map) {
            if (difference1.find(it.first) != difference1.end()) {
              kmer_maps[it.first] = "original";
            }
            if (difference2.find(it.first) != difference2.end()) {
              kmer_maps[it.first] = "mutated";
            }
          }
          return kmer_maps;
        }));
  }

  for (auto &thread_feature : thread_features) {
    auto thread_kmer_map = thread_feature.get();
    int count = 0;
    for (auto &it : thread_kmer_map) {
      std::cout << it.first << " " << it.second << std::endl;
      count++;
      if (count > 10) {
        break;
      }
    }
  }
  return kmer_map;
}

} // namespace kmer