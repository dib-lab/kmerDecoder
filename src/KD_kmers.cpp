#include "kmerDecoder.hpp"

void Kmers::seq_to_kmers(std::string& seq, std::vector<kmer_row>& kmers) {

    kmers.clear();
    kmers.reserve(seq.size());
    string filtered_seq;
    for (auto & c: seq) c = toupper(c);
    for (unsigned long i = 0; i < seq.size() - this->kSize + 1; i++) {
        kmer_row kmer;
        kmer.str = seq.substr(i, this->kSize);
        if (!valid_kmer(kmer.str)) continue;
        kmer.hash = this->hasher->hash(kmer.str);
        kmers.push_back(kmer);
    }

}

void Kmers::extractKmers() {

    bool SHORT_SEQ = false;

    for (int seqCounter = 0; seqCounter < Kmers::chunk_size && ((kseq_read(this->kseqObj)) >= 0); seqCounter++) {

        uint32_t seq_length = string(this->kseqObj->seq.s).size();

        if (seq_length < this->kSize) { SHORT_SEQ = true; continue; }

        std::string seq = kseqObj->seq.s;
        std::string id;
        id.append(kseqObj->name.s);
        if (kseqObj->comment.l) id.append(kseqObj->comment.s);

        this->kmers[id].reserve(seq.size());

        for (auto & c: seq) c = toupper(c);

        for (unsigned long i = 0; i < seq.size() - this->kSize + 1; i++) {
            kmer_row kmer;
            kmer.str = seq.substr(i, this->kSize);
            if (!valid_kmer(kmer.str)) continue;
            kmer.hash = this->hasher->hash(kmer.str);
            this->kmers[id].push_back(kmer);
        }

    }

    if ((unsigned int)this->kmers.size() != this->chunk_size && !SHORT_SEQ) {
        this->FILE_END = true;
    }

}

string Kmers::params_to_string() {
    string params;
    params = "kSize:" + to_string(this->kSize);
    return params;
}