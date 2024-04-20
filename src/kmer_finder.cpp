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

unordered_map<std::string, std::string>
parallel_finder(vector<unique_ptr<seq::Sequence>> &ref,
                unsigned int kmer_length, unsigned int num_threads,
                set<int> &difference1, set<int> &difference2) {

  unordered_map<std::string, std::string> kmer_map;
  vector<future<unordered_map<std::string, std::string>>> thread_features;
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
        thread_pool_->Submit([&ref, kmer_length, i, num_threads, part_size,
                              difference1, difference2]() {
          size_t start_idx = i * part_size;
          size_t end_idx = (i == num_threads - 1) ? ref.size()
                                                  : (start_idx + part_size);
          unordered_map<std::string, std::string> kmer_maps;
          vector<unique_ptr<seq::Sequence>> part_ref(
              make_move_iterator(ref.begin() + start_idx),
              make_move_iterator(ref.begin() + end_idx));
          
          for (auto &seq : part_ref) {
            auto set = kmer_maker_set(seq->genome, kmer_length);
            std::string name = seq->name;
            for (auto &it : set) {
              if (difference1.find(it) != difference1.end()) {
                kmer_maps[name] = "original";
                break;
              }
              else if (difference2.find(it) != difference2.end()) {
                kmer_maps[name] = "mutated";
                break;
              } else {
                kmer_maps[name] = "unknown";
              }
            }
          }

          return kmer_maps;
        }));
  }

  for (auto &thread_feature : thread_features) {
    auto thread_kmer_map = thread_feature.get();
    kmer_map.insert(thread_kmer_map.begin(), thread_kmer_map.end());
  }
  cout << endl << "Number of kmers: " << kmer_map.size() << endl;
  return kmer_map;
}

} // namespace kmer