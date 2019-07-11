#include "kmerDecoder.hpp"

void kmerDecoder::next_chunk(){
    seqan::clear(this->ids);
    seqan::clear(this->seqs);
    this->kmers.clear();
    seqan::readRecords(this->ids, this->seqs, *this->seqFileIn, this->chunk_size);
    this->seqan_end = seqan::atEnd(*this->seqFileIn);
    this->extractKmers();
}

flat_hash_map<std::string,std::vector<std::string>>* kmerDecoder::getKmers(){
    return &this->kmers;
}


bool kmerDecoder::end(){
    return this->seqan_end;
}