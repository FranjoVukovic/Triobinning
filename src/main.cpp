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

// Constants, adjust as needed
const int KMER_LENGTH = 21;
const int NUM_THREADS = 32;
const string paternalName = "PATERNAL";
const string maternalName = "MATERNAL";

int main() {
  // Path to the reference sequences, adjust as needed
  string pathRef = path + "paternal.fasta";
  string pathRef2 = path + "maternal.fasta";
  string pathRef3 = "/mnt/share1_Jabba/ftomas/ul_ec/ul_reads/subsampled/run2.output.filtered.subsampled.fasta";

  auto start = chrono::high_resolution_clock::now();

  // Read in the reference sequences and take out kmers
  auto p1 = bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(pathRef);
  auto ref = p1->Parse(-1);

  // Get all the kmers and write them to a file and a set
  auto map = kmer::parallel_kmer(ref, KMER_LENGTH, NUM_THREADS);
  ofstream outFile(path + "csv/outputPaternal.csv");
  set<int> keys1;
  for (auto &it : map) {
    keys1.insert(it.first);
    outFile << it.first << "," << it.second << endl;
  }
  map.clear();
  ref.clear();
  outFile.close();

  // Read in the reference sequences and take out kmers
  auto p2 = bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(pathRef2);
  auto ref2 = p2->Parse(-1);

  // Get all the kmers and write them to a file and a set
  auto map2 = kmer::parallel_kmer(ref2, KMER_LENGTH, NUM_THREADS);
  ofstream outFile2(path + "csv/outputMaternal.csv");
  set<int> keys2;
  for (auto &it : map2) {
    keys2.insert(it.first);
    outFile2 << it.first << "," << it.second << endl;
  }
  map2.clear();
  ref2.clear();
  outFile2.close();

  auto end = chrono::high_resolution_clock::now();

  // Calculate the time taken to read in the reference sequences and take out the unique kmers
  auto duration = chrono::duration_cast<chrono::seconds>(end - start).count();

  cout << "Time taken: " << duration << " seconds" << endl;
  cout << "Number of kmers: " << keys1.size() << endl;
  cout << "Number of kmers: " << keys2.size() << endl;

  // Reset the time
  auto currentTime = std::chrono::high_resolution_clock::now();
  auto currentTime_t = std::chrono::system_clock::to_time_t(currentTime);
  std::cout << "Start time: " << std::ctime(&currentTime_t) << std::endl;

  // Find the differences between the two sets of kmers; find the unique kmers
  set<int> diff;
  set_difference(keys1.begin(), keys1.end(), keys2.begin(), keys2.end(),
                 inserter(diff, diff.begin()));

  set<int> diff2;
  set_difference(keys2.begin(), keys2.end(), keys1.begin(), keys1.end(),
                 inserter(diff2, diff2.begin()));

  keys1.clear();
  keys2.clear();

  // See how much unique kmers there are for each parent
  cout << "Number of paternal unique kmers:" << diff.size() << endl << "Number of maternal unique kmers:" << diff2.size() << endl;

  // Read in the "child" sequences and see to which parent each kmer belongs
  auto p3 = bioparser::Parser<seq::Sequence>::Create<bioparser::FastaParser>(pathRef3);
  auto ref3 = p3->Parse(-1);

  auto myMap = kmer::parallel_finder(ref3, KMER_LENGTH, NUM_THREADS, diff, diff2);

  int fakePaternal = 0;
  int realPaternal = 0;
  int fakeMaternal = 0;
  int realMaternal = 0;
  int unknown = 0;

  for (auto &it : myMap) {
    // Get the name of the read
    auto name = it.first.substr(6, 8);
   
    // Check the real parent of the read
    if (name == paternalName) {
      if (it.second == 1) {
        realPaternal++;
      } else if (it.second == 0){
        unknown++;
      } else {
        fakePaternal++;
      }
    } else if (name == maternalName) {
      if (it.second == 2) {
        realMaternal++;
      } else if (it.second == 0){
        unknown++;
      } else if (it.second == 0){
        unknown++;
      } else {
        fakeMaternal++;
      }
    } 
  }

  // Print out the results
  cout << realPaternal << " " << fakePaternal << endl << fakeMaternal << " " << realMaternal << endl << unknown << endl;

  auto end2 = chrono::high_resolution_clock::now();
  auto duration2 = chrono::duration_cast<chrono::seconds>(end2 - end).count();

  cout << "Time taken to make differences: " << duration2 << " seconds" << endl;
  return 0;
}
