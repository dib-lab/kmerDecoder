#include "kmerDecoder.hpp"

void Kmers::seq_to_kmers(std::string &seq, std::vector<kmer_row> &kmers) {

    kmers.clear();

    kmers.reserve(seq.size());

    for (unsigned long i = 0; i < seq.size() - this->kSize + 1; i++) {
        kmer_row kmer;
        kmer.str = seq.substr(i, this->kSize);
        kmer.hash = this->hasher->hash(kmer.str);
        kmers.push_back(kmer);
    }

}

void Kmers::extractKmers() {

    for (int seqCounter = 0; seqCounter < this->chunk_size && ((kseq_read(this->kseqObj)) >= 0); seqCounter++) {

        uint32_t seq_length = string(this->kseqObj->seq.s).size();

        if (seq_length < this->kSize) continue;

        std::string seq = kseqObj->seq.s;
        std::string id = kseqObj->name.s;

        this->kmers[id].reserve(seq.size());

        for (unsigned long i = 0; i < seq.size() - this->kSize + 1; i++) {
            kmer_row kmer;
            kmer.str = seq.substr(i, this->kSize);
            kmer.hash = this->hasher->hash(kmer.str);
            this->kmers[id].push_back(kmer);
        }

    }

    if ((unsigned int) this->kmers.size() != this->chunk_size) {
        this->FILE_END = true;
    }

}

