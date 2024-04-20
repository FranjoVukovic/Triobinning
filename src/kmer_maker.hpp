#pragma once

#include "sequence.hpp"
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <set>

namespace kmer {

using namespace std;

set<unsigned int> kmer_maker_set(string genome, unsigned int kmer_length);

unordered_map<unsigned int, unsigned int>
kmer_maker(vector<unique_ptr<seq::Sequence>> &ref, unsigned int kmer_length);

unordered_map<unsigned int, unsigned int>
parallel_kmer(vector<unique_ptr<seq::Sequence>> &ref, unsigned int kmer_length,
              unsigned int num_threads);

} // namespace kmer