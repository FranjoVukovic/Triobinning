#include "kmer_finder.hpp"
#include "kmer_maker.hpp"
#include <fstream>
#include <future>
#include <iostream>
#include <set>
#include <string>
#include <thread_pool/thread_pool.hpp>
#include <unordered_map>
#include <vector>

namespace kmer {

static unsigned int THRESHOLD = 10;

unordered_map<std::string, unsigned int>
parallel_finder(vector<unique_ptr<seq::Sequence>> &ref,
                unsigned int kmer_length, unsigned int num_threads,
                set<int> &difference1, set<int> &difference2, unsigned int MIN, unsigned int MID) {

  // Create a thread pool
  vector<future<unordered_map<std::string, unsigned int>>> thread_features;
  shared_ptr<thread_pool::ThreadPool> thread_pool_ =
      make_shared<thread_pool::ThreadPool>(num_threads);

  const size_t part_size = ref.size() / num_threads;

  // Make num_threads threads and 
  for (uint32_t i = 0; i < num_threads; ++i) {
    thread_features.emplace_back(
        thread_pool_->Submit([&ref, kmer_length, i, num_threads, part_size,
                              &difference1, &difference2, MIN, MID]() {

          // Get the part of the reference                      
          size_t start_idx = i * part_size;
          size_t end_idx = (i == num_threads - 1) ? ref.size() : (start_idx + part_size);
          vector<unique_ptr<seq::Sequence>> part_ref(
              make_move_iterator(ref.begin() + start_idx),
              make_move_iterator(ref.begin() + end_idx));

          // Initialize the map
          unordered_map<std::string, unsigned int> kmer_maps;
          
          for (auto &seq : part_ref) {
            // Extract kmers from the read
            auto set = kmer_maker_set(seq->genome, kmer_length);
            
            std::string name = seq->name;
            int paternal = 0;
            int maternal = 0;

            // Compare the kmers with the parental kmers
            for (auto &it : set) {
              if (difference1.find(it) != difference1.end()) {
                paternal++;
              } else if (difference2.find(it) != difference2.end()) {
                maternal++;
              } 
            }

            // Assign the read to the parent with the most kmers
            if (paternal == 0 && maternal == 0)
              kmer_maps[name] = 0;       // Assign to unknown
            else if ((paternal >= MID && maternal <= MIN) || (paternal > maternal && maternal == 0) || (paternal/maternal > THRESHOLD))
              kmer_maps[name] = 1;      // Assign to paternal
            else if ((maternal >= MID && paternal <= MIN) || (maternal > paternal && paternal == 0) || (maternal/paternal > THRESHOLD))
              kmer_maps[name] = 2;      // Assign to maternal
            else
              kmer_maps[name] = 3;      // Assign to ambiguous
          }
          return kmer_maps;
        }));
  }

  unordered_map<std::string, unsigned int> kmer_map;

  // Combine the maps from the threads
  for (auto &thread_feature : thread_features) {
    auto thread_kmer_map = thread_feature.get();
    kmer_map.insert(thread_kmer_map.begin(), thread_kmer_map.end());
    thread_kmer_map.clear();
  }
  
  return kmer_map;
}

} // namespace kmer