#include "kmerDecoder.hpp"


void Kmers::extractKmers(){

    std::list<std::string> extracted_kmers;
    std::string id;
    std::string seq;

    for(unsigned seq_num = 0; seq_num < seqan::length(this->ids); seq_num++){

        seq = std::string((char*)seqan::toCString(this->seqs[seq_num]));
        id =  std::string((char*)seqan::toCString(this->ids[seq_num]));

        for (unsigned long i = 0; i < seq.size() - this->kSize + 1; i++)
        {
            // extracted_kmers.push_back(seq.substr(i, this->kSize));
            this->kmers[id].push_back(seq.substr(i, this->kSize));
        }

//        this->kmers[id] = extracted_kmers;
//        extracted_kmers.clear();
    }


}