#ifndef AAHASHER_HPP
#define AAHASHER_HPP

#include "hashUtils/hashutil.hpp"
#include <unordered_map>


/* ---------------------
 Class aaHasher_default: inherits from class:Hasher
 Hashing direction: Reversible
 Description:
    Used in hashing PROTEIN kmers with kSize <= 12.
    This is not a real hashing, it's a conversion from amino acids to {0:19} integers.
    Returns the non canonical 5-bit representation.
--------------------- */

class aaHasher_default : public Hasher {
protected:
    uint64_t kSize;
    std::unordered_map<char, uint8_t> complete_amino_to_int;
    std::unordered_map<uint8_t, char> complete_int_to_amino;

public:
    explicit aaHasher_default(uint64_t kSize);

    Hasher *clone() override { return new aaHasher_default(kSize); }

    uint64_t hash(const string &key) override;

    uint64_t hash(uint64_t key) override;

    string Ihash(uint64_t key) override;

    ~aaHasher_default(){}
};

class aaHasher_dayhoff : public aaHasher_default{
public:
    explicit aaHasher_dayhoff(uint64_t kSize);
    uint64_t hash(const string &key) override;
    string Ihash(uint64_t key) override;
    ~aaHasher_dayhoff(){}
};

#endif