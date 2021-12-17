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
*
*         Author:  Prashant Pandey (ppandey@cs.stonybrook.edu)
*                  Rob Patro (rob.patro@cs.stonybrook.edu)
*                  Rob Johnson (rob@cs.stonybrook.edu)
*   Organization:  Stony Brook University
*   Edited by: Mohamed Abuelanin (mabuelanin@gmail.com) UC Davis
*
* =====================================================================================
*/

#ifndef _HASHUTIL_H_
#define _HASHUTIL_H_

#include <sys/types.h>
#include <string>
#include <cstdlib>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include "Utils/kmer.h"
#include <algorithm>
#include <iostream>


using namespace std;

/* ---------------------
Class Hasher: Parent class
--------------------- */

class Hasher {
public:
    virtual uint64_t hash(const string &key) { return 0; };

    virtual uint64_t hash(uint64_t key) { return 0; };

    virtual string Ihash(uint64_t key) {
        throw logic_error("Reverese Hash function is not/ cannot be implemented for this hash function.");
    }

    virtual Hasher *clone() { return this; };

    virtual ~Hasher(){}
};


/* ---------------------
 Class bigKmerHasher: inherits from class:Hasher
 Hashing direction: Irreversible
 Description:
    Used in hashing large kmers with kSize > 32. It uses std::hash.
--------------------- */

class bigKmerHasher : public Hasher {
private:
    uint64_t kSize;
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

    ~bigKmerHasher(){};
};

/* ---------------------
 Class MumurHasher: inherits from class:Hasher
 Hashing direction: Irreversible
 Description:
    Used in hashing kmers with kSize <= 32. Supports non-ACGT kmers.
--------------------- */

class MumurHasher : public Hasher {
    using Hasher::hash;
private:
    uint64_t seed;
public:
    explicit MumurHasher(uint64_t Iseed) { seed = Iseed; }

    Hasher *clone() override { return new MumurHasher(seed); }

    uint64_t hash(const string & Skey) override;
    // uint64_t hash(uint64_t key) override;

    ~MumurHasher(){}
};

/* ---------------------
 Class IntegerHasher: inherits from class:Hasher
 Hashing direction: Reversible
 Description:
    Used in hashing kmers with kSize <= 32.
--------------------- */

class IntegerHasher : public Hasher {
private:
    uint64_t mask;
    uint64_t kSize;
public:
    explicit IntegerHasher(uint64_t kSize);

    Hasher *clone() override { return new IntegerHasher(kSize); }

    uint64_t hash(const string &key) override;

    uint64_t hash(uint64_t key) override;

    string Ihash(uint64_t key) override;

    ~IntegerHasher(){}
};

/* ---------------------
 Class noncanonical_IntegerHasher: inherits from class:IntegerHasher
 Hashing direction: Reversible
 Description:
    Used in hashing kmers with kSize <= 32.
    Returns the original non-canonical integerHash.
--------------------- */

class noncanonical_IntegerHasher : public Hasher {
private:
    uint64_t kSize;
    uint64_t mask;
  public:
    explicit noncanonical_IntegerHasher(uint64_t kSize){
        this->kSize = kSize;
        this->mask = BITMASK(2 * kSize);
    }

    Hasher *clone() override { return new noncanonical_IntegerHasher(kSize); }

    uint64_t hash(const string &key) override;
    uint64_t hash(uint64_t key) override;
    string Ihash(uint64_t key) override;
    
    ~noncanonical_IntegerHasher(){}
};


/* ---------------------
 Class TwoBitsHasher: inherits from class:Hasher
 Hashing direction: Reversible
 Description:
    Used in hashing kmers with kSize <= 32.
    This is not a real hashing, it's a conversion from A-C-G-T to 00-01-10-11.
    Returns the canonical two bit representation.
--------------------- */

class TwoBitsHasher : public Hasher {
private:
    uint64_t kSize;
public:
    explicit TwoBitsHasher(uint64_t kSize);

    Hasher *clone() override { return new TwoBitsHasher(kSize); }

    uint64_t hash(const string &key) override;

    uint64_t hash(uint64_t key) override;

    string Ihash(uint64_t key) override;

    ~TwoBitsHasher(){}
};

/* ---------------------
 Class noncanonical_TwoBitsHasher: inherits from class:TwoBitsHasher
 Hashing direction: Reversible
 Description:
    Used in hashing kmers with kSize <= 32.
    This is not a real hashing, it's a conversion from A-C-G-T to 00-01-10-11.
    Returns the original non-canonical two bit representation.
--------------------- */

class noncanonical_TwoBitsHasher : public TwoBitsHasher {
private:
    uint64_t kSize;
public:

    explicit noncanonical_TwoBitsHasher(uint64_t kSize) : TwoBitsHasher(kSize) {
        this->kSize = kSize;
    }

    Hasher *clone() override { return new noncanonical_TwoBitsHasher(kSize); }

    uint64_t hash(const string &key) override;

    ~noncanonical_TwoBitsHasher(){}
};

/* ---------------------
 Class QHasher: inherits from class:Hasher
 Hashing direction: Irreversible
 Description:
    Used in hashing kmers with kSize <= 32.
    Originally implemented to compress the index of kProcessor.
    Currently unused and could be removed in future releases.
--------------------- */

class QHasher : public Hasher {
private:
    uint64_t mask;
    uint64_t kSize;
    unsigned int Q = 28;
    unsigned int key_remainder_bits;

public:
    explicit QHasher(uint64_t kSize);

    QHasher(uint64_t kSize, int Q);

    Hasher *clone() override { return new QHasher(kSize, Q); }

    // To set the Q if not initialized
    void set_Q(int Q);

    uint64_t merge_Q_R(uint64_t &Q, uint64_t &R);

    void split_Q_R(uint64_t key, uint64_t &Q, uint64_t &R);

    uint64_t normal_hash(const string &key);

    uint64_t normal_hash(uint64_t key);

    uint64_t normal_Ihash(uint64_t key);


    uint64_t hash(const string &key) override;

    uint64_t hash(uint64_t key) override;

    string Ihash(uint64_t key) override;

    ~QHasher(){}
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

    uint64_t hash(const string &key) override {
        return fn(kmer::str_to_canonical_int(key));
    }

    uint64_t hash(uint64_t key) override {
        string kmerStr = kmer::int_to_str(key, kSize);
        return fn(key);
    }

    ~wrapperHasher(){}
};


#endif  // #ifndef _HASHUTIL_H_
