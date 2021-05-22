#include <doctest/doctest.h>
#include "kmerDecoder.hpp"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>
#include <string>

char letters[20] = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                    'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};

uint64_t get_kmer_int(string& kmer) {
  int kSize = kmer.size();
  Hasher* hasher = kmerDecoder::initHasher(protein_hasher, kSize);
  return hasher->hash(kmer);
}

string get_kmer_str(uint64_t& kmer, int& kSize) {
  Hasher* hasher = kmerDecoder::initHasher(protein_hasher, kSize);
  return hasher->Ihash(kmer);
}

TEST_CASE("Test aa-hashing") {
  for (int kSize = 1; kSize < 12; kSize++) {

    random_device rd;
    // Initialize Mersenne Twister pseudo-random number generator
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, 19);

    // Test x1000 times for this kmer    
    for (int j = 0; j < 1000; j++) {
      string aaKmer = "";
      // Generate random kmer
      for (int i = 0; i < kSize; i++) {
        aaKmer += letters[dis(gen)];
      }

      uint64_t aakmerHash = get_kmer_int(aaKmer);
      string calculated_aaKmer = get_kmer_str(aakmerHash, kSize);
      CHECK(aaKmer == calculated_aaKmer);
    }
  }
}