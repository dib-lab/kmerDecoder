#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <seqan/seq_io.h>
#include <parallel_hashmap/phmap.h>
#include <stdint.h>
#include "HashUtils/hashutil.hpp"

using phmap::flat_hash_map;


struct kmer_row {
    std::string str;
    uint64_t hash;
};

/* 
--------------------------------------------------------
                        InputModule:Parent
--------------------------------------------------------
*/

class kmerDecoder {

protected:
    unsigned int chunk_size;
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::CharString> seqs;
    flat_hash_map<std::string, std::vector<kmer_row>> kmers;
    std::string fileName;
    seqan::SeqFileIn seqFileIn;

    void initialize_seqan();

    bool seqan_end = false;

    virtual void extractKmers() = 0;


    // Mode 0: Murmar Hashing | Irreversible
    // Mode 1: Integer Hashing | Reversible | Full Hashing
    // Mode 2: TwoBitsHashing | Not considered hashing, just store the two bits representation


public:

    flat_hash_map<std::string, std::vector<kmer_row>> *getKmers();

    Hasher * hasher;

    int hash_mode = 0;
    bool canonical = true;

    virtual void seq_to_kmers(std::string &seq, std::vector<kmer_row> &kmers) = 0;

    virtual int get_kSize() = 0;

    bool end();

    void next_chunk();

    std::string get_filename();

    virtual void setHashingMode(int hash_mode, bool canonical = true) = 0;

    // hash single kmer
    uint64_t hash_kmer(std::string kmer_str) {
        return this->hasher->hash(kmer_str);
    }


    // Inverse hash single kmer
    std::string ihash_kmer(uint64_t kmer_hash) {
        return this->hasher->Ihash(kmer_hash);
    }

    static kmerDecoder *initialize_hasher(int kmer_size, int hash_mode = 1);


    static Hasher *create_hasher(int kmer_size, int hash_mode = 1) {

        switch (hash_mode) {
            case 0:
                return (new MumurHasher(2038074761));
            case 1:
                return (new IntegerHasher(kmer_size));
            case 2:
                return (new TwoBitsHasher(kmer_size));
            default:
                std::cerr << "Hashing mode : " << hash_mode << ", is not supported \n";
                std::cerr << "Mode 0: Murmar Hashing | Irreversible\n"
                             "Mode 1: Integer Hashing | Reversible\n"
                             "Mode 2: TwoBitsHashing | Not considered hashing, just store the two bits representation\n"
                          <<
                          "Default: Integer Hashing" << std::endl;
                exit(1);
        }

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

    void extractKmers();

public:

    explicit Kmers(int k_size, int hash_mode = 1) : kSize(k_size) {
        this->hasher = new IntegerHasher(kSize);
        this->hash_mode = 1;
        this->canonical = true;
        if (hash_mode != 1) {
            this->setHashingMode(hash_mode);
        }
    };

    Kmers(const std::string &filename, unsigned int chunk_size, int kSize) {
        this->kSize = kSize;
        this->fileName = filename;
        this->chunk_size = chunk_size;
        this->initialize_seqan();
        this->hasher = new IntegerHasher(kSize);
        this->hash_mode = 1;
        this->canonical = true;
    }

    void setHashingMode(int hash_mode, bool canonical = true) {

        this->hash_mode = hash_mode;
        this->canonical = canonical;

        if (hash_mode == 0) hasher = (new MumurHasher(2038074761));
        else if (hash_mode == 1) {
            if (canonical) hasher = (new IntegerHasher(kSize));
            else hasher = (new noncanonical_IntegerHasher(kSize));
        } else if (hash_mode == 2) {
            if (canonical) {
                hasher = (new TwoBitsHasher(kSize));
            } else {
                hasher = (new noncanonical_TwoBitsHasher(kSize));
            }
        } else {
            hasher = (new IntegerHasher(kSize));
        }

    }


    void seq_to_kmers(std::string &seq, std::vector<kmer_row> &kmers);


    int get_kSize() {
        return this->kSize;
    }
};


/* 
--------------------------------------------------------
                        Skipmers
--------------------------------------------------------
*/

class Skipmers : public kmerDecoder {
private:
    int m, n, k;
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
        this->hasher = new IntegerHasher(k);
        this->hash_mode = 1;
        this->canonical = true;
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
        this->fileName = filename;
        this->chunk_size = chunk_size;
        this->initialize_seqan();
        this->hasher = new IntegerHasher((int) k);
        this->hash_mode = 1;
        this->canonical = true;
    }

    void setHashingMode(int hash_mode, bool canonical = true) {
        this->hash_mode = hash_mode;
        this->canonical = canonical;
        if (hash_mode == 0) hasher = (new MumurHasher(2038074761));
        else if (hash_mode == 1) {
            if (canonical) hasher = (new IntegerHasher(k));
            else hasher = (new noncanonical_IntegerHasher(k));
        } else if (hash_mode == 2) {
            if (canonical) {
                hasher = (new TwoBitsHasher(k));
            } else {
                hasher = (new noncanonical_TwoBitsHasher(k));
            }
        } else {
            hasher = (new IntegerHasher(k));
        }

    }

    void seq_to_kmers(std::string &seq, std::vector<kmer_row> &kmers);

    int get_kSize() {
        return this->k;
    }

    virtual ~Skipmers() {}
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
        this->initialize_seqan();
        this->hasher = new IntegerHasher(k);
        this->hash_mode = 1;
        this->canonical = true;
    }

    Minimizers(int k, int w) {
        this->k = k;
        this->w = w;
        this->hasher = new IntegerHasher(k);
        this->hash_mode = 1;
        this->canonical = true;
    }

    void setHashingMode(int hash_mode, bool canonical = true) {
        this->hash_mode = hash_mode;
        this->canonical = canonical;
        if (hash_mode == 0) hasher = (new MumurHasher(2038074761));
        else if (hash_mode == 1) {
            if (canonical) hasher = (new IntegerHasher(k));
            else hasher = (new noncanonical_IntegerHasher(k));
        } else if (hash_mode == 2) {
            if (canonical) {
                hasher = (new TwoBitsHasher(k));
            } else {
                hasher = (new noncanonical_TwoBitsHasher(k));
            }
        } else {
            hasher = (new IntegerHasher(k));
        }
    }

    std::vector<mkmh_minimizer> getMinimizers(std::string &seq);


    void seq_to_kmers(std::string &seq, std::vector<kmer_row> &kmers);

    int get_kSize() {
        return this->k;
    }

    virtual ~Minimizers() {};

    static kmerDecoder *initialize_hasher(int kmer_size, int hash_mode = 1);
};