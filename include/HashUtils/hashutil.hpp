/*
* =====================================================================================
*
*       Filename:  hashutil.h
*
*    Description:
*
*        Version:  1.0
*        Created:  04/18/2016 04:49:32 PM
*       Revision:  none
*       Compiler:  gcc
*
*         Author:  Prashant Pandey (ppandey@cs.stonybrook.edu)
*                  Rob Patro (rob.patro@cs.stonybrook.edu)
*                  Rob Johnson (rob@cs.stonybrook.edu)
*   Organization:  Stony Brook University
*
* =====================================================================================
*/

#ifndef _HASHUTIL_H_
#define _HASHUTIL_H_

#include <sys/types.h>
#include <string>
#include <stdlib.h>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include "Utils/kmer.h"
#include <algorithm>
#include <boost/assign.hpp>
using namespace boost::assign;


using namespace std;

class Hasher {
public:
    virtual uint64_t hash(const string &key) { return 0; };

    virtual uint64_t hash(uint64_t key) { return 0; };

    virtual string Ihash(uint64_t key) {
        throw logic_error("Reverese Hash function is not/ cannot be implemented for this hash function.");
    }

    virtual Hasher *clone() { return this; };
};


class bigKmerHasher : public Hasher {
private:
    uint64_t kSize;
    map<char, char> m = map_list_of ('A', 'T') ('T', 'A') ('G', 'C') ('C', 'G');
    std::hash<string> hasher;
public:
    explicit bigKmerHasher(uint64_t kSize);

    string get_canonical_kmer(const string &kmer);

    uint64_t hash(const string &key) override;

    uint64_t hash(uint64_t key) override;

    string Ihash(uint64_t key) override {
        throw logic_error("Reverse Hash function is not/ cannot be implemented for this hash function.");
    }

    string Ihash(const string &key) {
        throw logic_error("Reverse Hash function is not/ cannot be implemented for this hash function.");
    }
};

class MumurHasher : public Hasher {
    using Hasher::hash;
private:
    uint64_t seed;
public:
    explicit MumurHasher(uint64_t Iseed) { seed = Iseed; }

    Hasher *clone() override { return new MumurHasher(seed); }

    uint64_t hash(string kmer);
};

class IntegerHasher : public Hasher {
private:
    uint64_t mask;
    uint64_t kSize;
public:
    IntegerHasher(uint64_t kSize);

    Hasher *clone() override { return new IntegerHasher(kSize); }

    uint64_t hash(const string &key) override;

    uint64_t hash(uint64_t key) override;

    string Ihash(uint64_t key) override;
};

class QHasher : public Hasher {
private:
    uint64_t mask;
    uint64_t kSize;
    unsigned int Q = 28;
    unsigned int key_remainder_bits;

public:
    explicit QHasher(uint64_t kSize);

    QHasher(uint64_t kSize, int _Q);

    Hasher *clone() override { return new QHasher(kSize, Q); }

    // To set the Q if not initialized
    void set_Q(int _Q);

    uint64_t merge_Q_R(uint64_t &_Q, uint64_t &R);

    void split_Q_R(uint64_t key, uint64_t &_Q, uint64_t &R);

    uint64_t normal_hash(const string &key);

    uint64_t normal_hash(uint64_t key);

    uint64_t normal_Ihash(uint64_t key);


    uint64_t hash(const string &key) override;

    uint64_t hash(uint64_t key) override;

    string Ihash(uint64_t key) override;
};

// TwoBitsHasher Canonical

class TwoBitsHasher : public Hasher {
private:
    uint64_t kSize;
public:
    explicit TwoBitsHasher(uint64_t kSize);

    Hasher *clone() override { return new TwoBitsHasher(kSize); }

    uint64_t hash(const string &key) override;

    uint64_t hash(uint64_t key) override;

    string Ihash(uint64_t key) override;
};


// TwoBitsHasher nonCanonical

class noncanonical_TwoBitsHasher : public TwoBitsHasher {
private:
    uint64_t kSize;
public:

    noncanonical_TwoBitsHasher(uint64_t kSize) : TwoBitsHasher(kSize) {
        this->kSize = kSize;
    }

    Hasher *clone() override { return new noncanonical_TwoBitsHasher(kSize); }

    uint64_t hash(const string &key) override;
};


class noncanonical_IntegerHasher : public IntegerHasher {
private:
    uint64_t kSize;
    uint64_t mask;
public:
    explicit noncanonical_IntegerHasher(uint64_t kSize) : IntegerHasher(kSize) {
        this->kSize = kSize;
        this->mask = BITMASK(2 * kSize);
    }

    Hasher *clone() override { return new noncanonical_IntegerHasher(kSize); }

    uint64_t hash(const string &key) override;
};


template<class hashFnType>
class wrapperHasher : public Hasher {
private:
    hashFnType fn;
    uint64_t kSize;
public:
    wrapperHasher(hashFnType fnn, uint64_t kSize)
            : fn(fnn), kSize(kSize) {}

    Hasher *clone() override {
        return new wrapperHasher(fn, kSize);
    }

    uint64_t hash(const string &key) {
        return fn(kmer::str_to_canonical_int(key));
    }

    uint64_t hash(uint64_t key) override {
        string kmerStr = kmer::int_to_str(key, kSize);
        return fn(key);
    }

};


#endif  // #ifndef _HASHUTIL_H_
