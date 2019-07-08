#include "kmerDecoder.hpp"


void Skipmers::extractKmers()
{
    std::string ORF_SEQ = "";

    std::vector<std::string> extracted_kmers;
    std::string id;
    std::string seq;

    for(int seq_num =0; seq_num < this->chunk_size; seq_num++){

        id = std::string((char*)seqan::toCString(this->ids[seq_num]));
        seq = std::string((char*)seqan::toCString(this->seqs[seq_num]));

        for(int start = 0; start < 3; start++){

            for(unsigned long i = start; i < seq.size(); i+=this->n){
                ORF_SEQ.append(seq.substr(i, this->m));
            }

            for(unsigned long j = 0; j < ORF_SEQ.size() - this->k + 1; j++){
                extracted_kmers.push_back(ORF_SEQ.substr(j,this->k));
            }

            ORF_SEQ.clear();
        }

        this->kmers[id] = extracted_kmers;
        extracted_kmers.clear();
    }

}