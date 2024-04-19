#pragma once

#include "sequence.hpp"
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>


namespace kmer {

unordered_map<unsigned int, std::string>
parallel_finder(vector<unique_ptr<seq::Sequence>> &ref,
                unsigned int kmer_length, unsigned int num_threads,
                set<int> &difference1, set<int> &difference2);

} // namespace kmer