#include "hashUtils/aaHasher.hpp"

aaHasher_default::aaHasher_default(uint64_t kSize) {
  this->kSize = kSize;
  this->complete_amino_to_int = {
      {'A', 0},
      {'R', 1},
      {'N', 2},
      {'D', 3},
      {'C', 4},
      {'Q', 5},
      {'E', 6},
      {'G', 7},
      {'H', 8},
      {'I', 9},
      {'L', 10},
      {'K', 12},
      {'M', 13},
      {'F', 14},
      {'P', 15},
      {'S', 16},
      {'T', 17},
      {'W', 18},
      {'Y', 19},
      {'V', 20}
  };

  this->complete_int_to_amino = {
      {0 ,'A'},
      {1 ,'R'},
      {2 ,'N'},
      {3 ,'D'},
      {4 ,'C'},
      {5 ,'Q'},
      {6 ,'E'},
      {7 ,'G'},
      {8 ,'H'},
      {9 ,'I'},
      {10 ,'L'},
      {12 ,'K'},
      {13 ,'M'},
      {14 ,'F'},
      {15 ,'P'},
      {16 ,'S'},
      {17 ,'T'},
      {18 ,'W'},
      {19 ,'Y'},
      {20 ,'V'},
  };
}

uint64_t aaHasher_default::hash(const string &aa_str_kmer) {

    uint64_t hash_val = 0;
    for (auto const &ch : aa_str_kmer) {
        uint8_t aa_int = complete_amino_to_int[ch];
        hash_val = hash_val | aa_int;
        hash_val = hash_val << 5;
    }

    return hash_val >> 5;
}

uint64_t aaHasher_default::hash(uint64_t key) {
    return key;
}

string aaHasher_default::Ihash(uint64_t aa_int_kmer) {

    string aa_str_kmer;
    
    for (auto i = (uint8_t)this->kSize; i > 0; i--) {
        uint8_t base = (aa_int_kmer >> (i * 5 - 5)) & 31ULL;
        aa_str_kmer.push_back(complete_int_to_amino[base]);
    }

    return aa_str_kmer;
}


// DAYHOFF

aaHasher_dayhoff::aaHasher_dayhoff(uint64_t kSize) : aaHasher_default(kSize) {
    
  this->kSize = kSize;
  this->complete_amino_to_int = {
      {'C', 0},
      {'G', 1}, {'S', 1}, {'T', 1}, {'A', 1}, {'P', 1},
      {'D', 2}, {'E', 2}, {'N', 2}, {'Q', 2},
      {'R', 3}, {'H', 3},{'K', 3},
      {'L', 4}, {'V', 4}, {'M', 4}, {'I', 4},
      {'Y', 5}, {'F', 5}, {'W', 5},
  };

  this->complete_int_to_amino = {};


}

uint64_t aaHasher_dayhoff::hash(const string &aa_str_kmer) {

    uint64_t hash_val = 0;
    for (auto const &ch : aa_str_kmer) {
        uint8_t aa_int = complete_amino_to_int[ch];
        hash_val = hash_val | aa_int;
        hash_val = hash_val << 3;
    }

    return hash_val >> 3;
}

string aaHasher_dayhoff::Ihash(uint64_t key)
{
    throw "Can't reverse hash Dayhoff encoding";
}
