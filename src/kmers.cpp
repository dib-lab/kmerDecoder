#include "inputModule.hpp"


std::list<std::string>* Kmers::getKmers(std::string &seq){

    this->kmers.clear();

    for (int i = 0; i < seq.size() - this->kSize + 1; i++)
    {
        this->kmers.push_back(seq.substr(i, this->kSize));
    }

    return &this->kmers;

}