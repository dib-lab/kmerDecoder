#include "inputModule.hpp"


std::list<std::string>* Skipmers::getKmers(std::string &seq)
{   
    this->kmers.clear();
    std::string ORF_SEQ = "";

    for(int start = 0; start < 3; start++){

        for(int i = start; i < seq.size(); i+=this->n){
            ORF_SEQ.append(seq.substr(i, this->m));
        }

        for(int j = 0; j < ORF_SEQ.size() - this->k + 1; j++){
            this->kmers.push_back(ORF_SEQ.substr(j,this->k));
        }

        ORF_SEQ.clear();
    }

    return &this->kmers;

}