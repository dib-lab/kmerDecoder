#include "kmerDecoder.hpp"

void kmerDecoder::next_chunk(){
    seqan::clear(this->ids);
    seqan::clear(this->seqs);
    this->kmers.clear();
    seqan::readRecords(this->ids, this->seqs, this->seqFileIn, this->chunk_size);
    this->seqan_end = seqan::atEnd(this->seqFileIn);
    this->extractKmers();
}

flat_hash_map<std::string,std::vector<kmer_row>>* kmerDecoder::getKmers(){
    return &this->kmers;
}

void kmerDecoder::initialize_seqan(){

    if (!seqan::open(this->seqFileIn, seqan::toCString(this->fileName)))
    {
        std::cerr << "ERROR: Could not open the file.\n";
        exit(1);
    }

}

bool kmerDecoder::end(){
    return this->seqan_end;
}

std::string kmerDecoder::get_filename(){
    return this->fileName;
}


kmerDecoder * kmerDecoder::initialize_hasher(int kmer_size, int hash_mode){
    return new Kmers(kmer_size, hash_mode);
}