#include "kmerDecoder.hpp"


void Kmers::extractKmers(){

    std::string id;
    std::string seq;

    for(unsigned seq_num = 0; seq_num < seqan::length(this->ids); seq_num++){
        seq = std::string((char*)seqan::toCString(this->seqs[seq_num]));
        id =  std::string((char*)seqan::toCString(this->ids[seq_num]));

        this->kmers[id].reserve(seq.size());

        for (unsigned long i = 0; i < seq.size() - this->kSize + 1; i++)
        {
            this->kmers[id].push_back(seq.substr(i, this->kSize));
        }

    }


}