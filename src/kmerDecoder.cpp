#include "kmerDecoder.hpp"

void kmerDecoder::next_chunk(){
    this->kmers.clear();
    this->extractKmers();
}

flat_hash_map<std::string,std::vector<kmer_row>>* kmerDecoder::getKmers(){
    return &this->kmers;
}


void kmerDecoder::initialize_kSeq(){
    fp = gzopen(this->fileName.c_str(), "r");
    kseqObj = kseq_init(fp);
}

bool kmerDecoder::end(){
    return this->FILE_END;
}

std::string kmerDecoder::get_filename(){
    return this->fileName;
}


kmerDecoder * kmerDecoder::initialize_hasher(int kmer_size, int hash_mode){
    return new Kmers(kmer_size, hash_mode);
}