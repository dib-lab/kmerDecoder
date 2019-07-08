#include "kmerDecoder.hpp"

void kmerDecoder::next_chunk(){
    seqan::clear(this->ids);
    seqan::clear(this->seqs);
    this->kmers.clear();
    seqan::readRecords(this->ids, this->seqs, this->seqFileIn, this->chunk_size);
    this->seqan_end = seqan::atEnd(this->seqFileIn);
    this->extractKmers();
}


void kmerDecoder::initialize_seqan(){

    if (!seqan::open(this->seqFileIn, seqan::toCString(this->seqFileName)))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        exit(1);
    }

    this->next_chunk();
}

void kmerDecoder::set_fileName(std::string filename) {
    this->seqFileName = (seqan::CharString)filename;
}

flat_hash_map<std::string,std::vector<std::string>>* kmerDecoder::getKmers(){
    return &this->kmers;
}

kmerDecoder &kmerDecoder::operator++(){
    this->next_chunk();
}

bool kmerDecoder::end(){
    return this->seqan_end;
}