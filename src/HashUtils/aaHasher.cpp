//
// Created by mabuelanin on 2/17/21.
//

#include "HashUtils/aaHasher.hpp"


std::unordered_map<char, uint8_t> complete_amino_to_int{
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

std::unordered_map<uint8_t, char> complete_int_to_amino{
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

aaHasher::aaHasher(uint64_t kSize) {
    this->kSize = kSize;
}

uint64_t aaHasher::hash(const string &aa_str_kmer) {

    uint64_t hash_val = 0;
    for (auto const &ch : aa_str_kmer) {
        uint8_t aa_int = complete_amino_to_int[ch];
        hash_val = hash_val | aa_int;
        hash_val = hash_val << 5;
    }

    return hash_val >> 5;
}

uint64_t aaHasher::hash(uint64_t key) {
    return key;
}

string aaHasher::Ihash(uint64_t aa_int_kmer) {

    string aa_str_kmer;
    
    for (uint8_t i = (uint8_t)this->kSize; i > 0; i--) {
        uint8_t base = (aa_int_kmer >> (i * 5 - 5)) & 31ULL;
        aa_str_kmer.push_back(complete_int_to_amino[base]);
    }

    return aa_str_kmer;
}
