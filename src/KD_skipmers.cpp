#include "kmerDecoder.hpp"

void Skipmers::seq_to_kmers(std::string & seq, std::vector <kmerRow> & kmers){

    kmers.clear();
    kmers.reserve(seq.size());

    std::string ORF_SEQ = "";
    for(auto const & start : this->ORFs){

        for(unsigned long i = start; i < seq.size(); i+=this->n){
                ORF_SEQ.append(seq.substr(i, this->m));
        }

        for(unsigned long j = 0; j < ORF_SEQ.size() - this->k + 1; j++){
            kmerRow kmer;
            kmer.str = ORF_SEQ.substr(j,this->k);
            kmer.hash = this->hasher->hash(kmer.str);
            kmers.push_back(kmer);
        }
        
        ORF_SEQ.clear();
    }

}

void Skipmers::extractKmers()
{
    std::string ORF_SEQ = "";
    std::string id;
    std::string seq;

    for(unsigned long seq_num =0; seq_num < seqan::length(this->ids); seq_num++){

        id = std::string((char*)seqan::toCString(this->ids[seq_num]));
        seq = std::string((char*)seqan::toCString(this->seqs[seq_num]));
        
        this->kmers[id].reserve(seq.size());

        for(auto const & start : this->ORFs){

            for(unsigned long i = start; i < seq.size(); i+=this->n){
                ORF_SEQ.append(seq.substr(i, this->m));
            }

            for(unsigned long j = 0; j < ORF_SEQ.size() - this->k + 1; j++){
                kmerRow kmer;
                kmer.str = ORF_SEQ.substr(j,this->k);
                kmer.hash = this->hasher->hash(kmer.str);
                this->kmers[id].push_back(kmer);
            }

            ORF_SEQ.clear();
        }
    }

}