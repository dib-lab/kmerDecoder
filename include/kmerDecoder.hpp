#ifndef KMERDECODER_HPP
#define KMERDECODER_HPP

#include <kseq/kseq.h>
#include <parallel_hashmap/phmap.h>
#include <zlib.h>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "hashUtils/aaHasher.hpp"

KSEQ_INIT(gzFile, gzread)

using phmap::flat_hash_map;

// using KMERS_MAP = parallel_flat_hash_map<std::string, std::vector<kmer_row>,
//         phmap::priv::hash_default_hash<std::string>,
//         phmap::priv::hash_default_eq<std::string>,
//         std::allocator<std::pair<const std::string,
//         std::vector<kmer_row>>>,
//         1
//         // ,std::mutex
//         >;

enum readingModes{
  KMERS = 1,
  SKIPMERS = 2,
  MINIMIZERS = 3,
  PROTEIN = 4,
};

enum hashingModes{
  mumur_hasher = 0,
  integer_hasher = 1,
  TwoBits_hasher = 2,
  nonCanonicalInteger_Hasher = 3,
  protein_hasher = 4,
  proteinDayhoff_hasher = 5,
};

struct kmer_row {
    std::string str;
    uint64_t hash;
};


bool valid_kmer(const string& kmer);


/* 
--------------------------------------------------------
                        InputModule:Parent
--------------------------------------------------------
*/



class kmerDecoder {

protected:
    unsigned int chunk_size{};

    flat_hash_map<std::string, std::vector<kmer_row>> kmers;
    std::string fileName;
    gzFile fp{};
    kseq_t *kseqObj{};


    void initialize_kSeq();

    bool FILE_END = false;

    virtual void extractKmers() = 0;


public:
    static Hasher* initHasher(hashingModes HM, int kSize);

    flat_hash_map<std::string, std::vector<kmer_row>> *getKmers();

    static kmerDecoder *getInstance(readingModes RM, hashingModes HM, map<string, int> params);
    static kmerDecoder * getInstance(string fileName, int chunkSize, readingModes RM, hashingModes HM, map<string, int> params);

    Hasher * hasher{};

    hashingModes hash_mode = integer_hasher;
    readingModes slicing_mode;

    virtual void seq_to_kmers(std::string &seq, std::vector<kmer_row> &kmers) = 0;

    virtual int get_kSize() = 0;

    bool end() const;

    void next_chunk();

    std::string get_filename();

    virtual void setHashingMode(hashingModes HM, int kSize) = 0;

    static std::map<std::string, int> string_to_params(std::string const& s);
    virtual string params_to_string() = 0;

    // hash single kmer
    uint64_t hash_kmer(const std::string & kmer_str) {
        return this->hasher->hash(kmer_str);
    }


    // Inverse hash single kmer
    std::string ihash_kmer(uint64_t kmer_hash) {
        return this->hasher->Ihash(kmer_hash);
    }

    virtual ~kmerDecoder(){
        kseq_destroy(this->kseqObj);
        gzclose(this->fp);
        this->kmers.clear();
    }

};


/* 
--------------------------------------------------------
                        Default Kmers
--------------------------------------------------------
*/


class Kmers : public kmerDecoder {

private:
    unsigned kSize{};
    void extractKmers() override;

public:

    explicit Kmers(int k_size, hashingModes HM = integer_hasher) : kSize(k_size) {
        this->hasher = kmerDecoder::initHasher(HM, kSize);
        this->slicing_mode = KMERS;
        this->hash_mode = HM;
    };

    Kmers(const std::string &filename, unsigned int chunk_size, int kSize, hashingModes HM = integer_hasher) {
        this->kSize = kSize;
        this->fileName = filename;
        this->chunk_size = chunk_size;
        this->initialize_kSeq();
        this->hasher = kmerDecoder::initHasher(HM, kSize);
        this->hash_mode = HM;
        this->slicing_mode = KMERS;
    }

    void setHashingMode(hashingModes HM, int _kSize = 0) {
        if(_kSize) this->kSize = _kSize;
        this->hash_mode = HM;
        delete this->hasher;
        this->hasher = kmerDecoder::initHasher(HM, this->kSize);
    }


    void seq_to_kmers(std::string &seq, std::vector<kmer_row> &kmers) override;
    string params_to_string() override;


    int get_kSize() {
        return this->kSize;
    }

    ~Kmers() {}

};


/* 
--------------------------------------------------------
                        Skipmers
--------------------------------------------------------
*/

class Skipmers : public kmerDecoder {
private:
    int m, n, k, S;
    std::vector<int> ORFs = {0, 1, 2};

    void extractKmers();

public:

    Skipmers(uint8_t m, uint8_t n, uint8_t k, int ORF = 0) {
        if (n < 1 || n < m || k < m || k % m != 0) {
            std::cerr << "Error: invalid skip-mer shape!"
                      << "Conditions: 0 < m <= n < k & k must be multiple of m" << std::endl;

            exit(1);
        }

        if (ORF) {
            this->ORFs.clear();
            this->ORFs.push_back(ORF - 1);
        }

        this->m = m;
        this->n = n;
        this->k = k;
        this->S = k;
        this->S = S + ((S - 1) / this->m) * (this->n - this->m);
        this->hasher = kmerDecoder::initHasher(integer_hasher, k);
        this->hash_mode = integer_hasher;
        this->slicing_mode = SKIPMERS;
    }

    Skipmers(const std::string &filename, unsigned int chunk_size, uint8_t m, uint8_t n, uint8_t k, int ORF = 0) {
        if (n < 1 or n < m || k < m || k % m != 0) {
            std::cerr << "Error: invalid skip-mer shape!"
                      << "Conditions: 0 < m <= n < k & k must be multiple of m" << std::endl;
            exit(1);
        }

        if (ORF) {
            this->ORFs.clear();
            this->ORFs.push_back(ORF - 1);
        }

        this->m = m;
        this->n = n;
        this->k = k;
        this->S = k;
        this->S = S + ((S - 1) / this->m) * (this->n - this->m);
        this->fileName = filename;
        this->chunk_size = chunk_size;
        this->initialize_kSeq();
        this->hasher = kmerDecoder::initHasher(integer_hasher, k);
        this->hash_mode = integer_hasher;
        this->slicing_mode = SKIPMERS;
    }

    void setHashingMode(hashingModes HM, int _kSize = 0) {
        if(_kSize) this->k = _kSize;
        this->hash_mode = HM;
        delete this->hasher;
        this->hasher = kmerDecoder::initHasher(HM, this->k);
    }
    
    void seq_to_kmers(std::string &seq, std::vector<kmer_row> &kmers);
    string params_to_string() override;

    int get_kSize() {
        return this->k;
    }

    ~Skipmers() {}
};


/* 
--------------------------------------------------------
                        Minimizers
--------------------------------------------------------
*/

typedef struct mkmh_minimizer {
    uint64_t pos;
    uint32_t length;
    std::string seq;

    bool operator<(const mkmh_minimizer &rhs) const { return seq < rhs.seq; };
} mkmh_minimizer;


class Minimizers : public kmerDecoder {
private:
    int k, w;

    void extractKmers();

    struct mkmh_kmer_list_t {
        char **kmers;
        int length;
        int k;

        mkmh_kmer_list_t() {

        };

        mkmh_kmer_list_t(int length, int k) {
            this->length = length;
            this->k = k;
            this->kmers = new char *[length];
        };

        ~mkmh_kmer_list_t() {
            for (int i = 0; i < this->length; ++i) {
                delete[] this->kmers[i];
            }
            delete[] this->kmers;
        };
    };

protected:
    std::vector<mkmh_minimizer> kmer_tuples(std::string seq, int k);

    mkmh_kmer_list_t kmerize(char *seq, int seq_len, int k);

    std::vector<std::string> kmerize(std::string seq, int k);

    void kmerize(char *seq, const int &seq_len, const int &k, char **kmers, int &kmer_num);

    template<typename T>
    std::vector<T> v_set(std::vector<T> kmers);

public:
    Minimizers(const std::string &filename, unsigned int chunk_size, int k, int w) {
        this->k = k;
        this->w = w;
        this->fileName = filename;
        this->chunk_size = chunk_size;
        this->initialize_kSeq();
        this->hasher = kmerDecoder::initHasher(integer_hasher, k);
        this->hash_mode = integer_hasher;
        this->slicing_mode = MINIMIZERS;
    }

    Minimizers(int k, int w) {
        this->k = k;
        this->w = w;
        this->hasher = kmerDecoder::initHasher(integer_hasher, k);
        this->hash_mode = integer_hasher;
        this->slicing_mode = MINIMIZERS;
    }

    void setHashingMode(hashingModes HM, int _kSize = 0) {
        if(_kSize) this->k = _kSize;
        this->hash_mode = HM;
        delete this->hasher;
        this->hasher = kmerDecoder::initHasher(HM, this->k);
    }

    std::vector<mkmh_minimizer> getMinimizers(std::string &seq);


    void seq_to_kmers(std::string &seq, std::vector<kmer_row> &kmers);
    string params_to_string() override;

    int get_kSize() {
        return this->k;
    }

    ~Minimizers(){
        kseq_destroy(this->kseqObj);
        gzclose(this->fp);
        this->kmers.clear();
    }

};


/* 
--------------------------------------------------------
                        AA Kmers (Protein Seqs)
--------------------------------------------------------
*/


class aaKmers : public kmerDecoder {

private:
    unsigned kSize{};

    void extractKmers() override;

public:

    explicit aaKmers(int k_size, hashingModes HM = protein_hasher) : kSize(k_size) {

        if((kSize > 11 && HM == protein_hasher) || (kSize > 19 && HM == proteinDayhoff_hasher)){
            throw "can't use aaKmer > 11";
        }

        this->hasher = kmerDecoder::initHasher(HM, kSize);
        this->slicing_mode = PROTEIN;
        this->hash_mode = HM;
    };

    aaKmers(const std::string &filename, unsigned int chunk_size, int kSize, hashingModes HM = protein_hasher) {

        if(kSize > 11){
            throw "can't use aaKmer > 11";
        }

        this->kSize = kSize;
        this->fileName = filename;
        this->chunk_size = chunk_size;
        this->initialize_kSeq();
        this->hasher = kmerDecoder::initHasher(HM, kSize);
        this->hash_mode = HM;
        this->slicing_mode = PROTEIN;
    }

    void seq_to_kmers(std::string &seq, std::vector<kmer_row> &kmers) override;
    string params_to_string() override;


    void setHashingMode(hashingModes HM, int _kSize = 0) {
        if(_kSize) this->kSize = _kSize;
        this->hash_mode = HM;
        delete this->hasher;
        this->hasher = kmerDecoder::initHasher(HM, this->kSize);
    }

    int get_kSize() {
        return this->kSize;
    }

    ~aaKmers() override{
        delete this->hasher;
        kseq_destroy(this->kseqObj);
        gzclose(this->fp);
        this->kmers.clear();
    }

};
#endif