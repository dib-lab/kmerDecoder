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
    KS_FULL_COMMENT = true; // Allow retrieving full header from kseq
}

bool kmerDecoder::end() const{
    return this->FILE_END;
}

std::string kmerDecoder::get_filename(){
    return this->fileName;
}

Hasher* kmerDecoder::initHasher(hashingModes HM, int kSize)
{
    switch (HM) {
        case integer_hasher:return new IntegerHasher(kSize);
        case mumur_hasher:return new MumurHasher(2038074761);
        case TwoBits_hasher:return new TwoBitsHasher(kSize);
        case nonCanonicalInteger_Hasher:return new noncanonical_IntegerHasher(kSize);
        case protein_hasher:return new aaHasher_default(kSize);
        case proteinDayhoff_hasher:return new aaHasher_dayhoff(kSize);
        default: throw "unknown hashing mode";
    }
}
