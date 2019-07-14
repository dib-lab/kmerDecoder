#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <seqan/seq_io.h>
#include <parallel_hashmap/phmap.h>

using phmap::flat_hash_map;

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
    flat_hash_map<std::string, std::vector<std::string>> kmers;
    std::string fileName;
    seqan::SeqFileIn seqFileIn;
    void initialize_seqan();
    bool seqan_end = false;
    virtual void extractKmers() = 0;

public:

    flat_hash_map<std::string, std::vector<std::string>> *getKmers();

    
    virtual void seq_to_kmers(std::string & seq, std::vector <std::string> & kmers) = 0;
    virtual int get_kSize() = 0;

    bool end();

    void next_chunk();

    std::string get_filename();

};


/* 
--------------------------------------------------------
                        Default Kmers
--------------------------------------------------------
*/


class Kmers : public kmerDecoder {

private:
    unsigned kSize;
    void extractKmers();

public:

    Kmers(const std::string & filename, unsigned int chunk_size, int kSize) {
        this->kSize = kSize;
        this->fileName = filename;
        this->chunk_size = chunk_size;
        this->initialize_seqan();
    }

    explicit Kmers(int k_size) : kSize(k_size) {};

    void seq_to_kmers(std::string & seq, std::vector <std::string> & kmers);


    int get_kSize(){
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
    void extractKmers();

public:

    Skipmers(uint8_t m, uint8_t n, uint8_t k){
        if (n < 1 || n < m || k < m || k % m != 0) {
            std::cerr << "Error: invalid skip-mer shape!"
                      << "Conditions: 0 < m <= n < k & k must be multiple of m" << std::endl;

            exit(1);
        }

        this->m = m;
        this->n = n;
        this->k = k;
    }

    Skipmers(const std::string & filename, unsigned int chunk_size, uint8_t m, uint8_t n, uint8_t k) {
        if (n < 1 or n < m || k < m || k % m != 0) {
            std::cerr << "Error: invalid skip-mer shape!"
                      << "Conditions: 0 < m <= n < k & k must be multiple of m" << std::endl;
            exit(1);
        }

        this->m = m;
        this->n = n;
        this->k = k;
        this->fileName = filename;
        this->chunk_size = chunk_size;
        this->initialize_seqan();
    }


    void seq_to_kmers(std::string & seq, std::vector <std::string> & kmers);

    int get_kSize(){
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
    Minimizers(const std::string & filename, unsigned int chunk_size, int k, int w) {
        this->k = k;
        this->w = w;
        this->fileName = filename;
        this->chunk_size = chunk_size;
        this->initialize_seqan();
    }

    Minimizers(int k, int w) {
        this->k = k;
        this->w = w;
    }

    std::vector<mkmh_minimizer> getMinimizers(std::string &seq);


    void seq_to_kmers(std::string & seq, std::vector <std::string> & kmers);

    int get_kSize(){
        return this->k;
    }

    virtual ~Minimizers() {};
};