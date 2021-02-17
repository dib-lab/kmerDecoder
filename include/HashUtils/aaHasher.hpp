//
// Created by mabuelanin on 2/17/21.
//

#ifndef KMERDECODER_AAHASHER_HPP
#define KMERDECODER_AAHASHER_HPP

#include "hashutil.hpp"
#include <unordered_map>


/* ---------------------
 Class aaHasher: inherits from class:Hasher
 Hashing direction: Reversible
 Description:
    Used in hashing PROTEIN kmers with kSize <= 12.
    This is not a real hashing, it's a conversion from amino acids to {0:19} integers.
    Returns the non canonical 5-bit representation.
--------------------- */

class aaHasher : public Hasher {
private:
    uint64_t kSize;

public:
    explicit aaHasher(uint64_t kSize);

    Hasher *clone() override { return new aaHasher(kSize); }

    uint64_t hash(const string &key) override;

    uint64_t hash(uint64_t key) override;

    string Ihash(uint64_t key) override;
};


#endif //KMERDECODER_AAHASHER_HPP
