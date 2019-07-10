#include <vector>
#include <queue>
#include <set>
#include <unordered_map>
#include <string>
#include <cstring>
#include <sstream>
#include <locale>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <omp.h>
#include <assert.h>
#include <bitset>
#include "kmerDecoder.hpp"

using std::vector;
using std::queue;
using std::string;
using std::set;
using std::unordered_map;

/* 
  --------------------------------------------------------
              Derived Class : Minimizers
  --------------------------------------------------------
*/


/* Returns the forward kmers of a sequence */
vector<string> Minimizers::kmerize(string seq, int k)
{
    vector<string> ret(seq.length() - k, "");

#pragma omp parallel for
    for (int i = 0; i < seq.length() - k; i++)
    {
        string s = seq.substr(i, k);
        //#pragma omp atomic read
        ret[i] = s;
        //ret.push_back(s);
        //ret.push_back(reverse(reverse_complement(s)));
    }
    return ret;
};

Minimizers::mkmh_kmer_list_t Minimizers::kmerize(char *seq, int seq_len, int k)
{
    mkmh_kmer_list_t ret;
    ret.kmers = new char *[seq_len - k];
    ret.k = k;
    ret.length = seq_len - k;

    for (int i = 0; i < ret.length; ++i)
    {
        char *km = new char[k + 1];
        memcpy(km, seq + i, k);
        ret.kmers[i] = new char[k + 1];
        ret.kmers[i] = km;
        ret.kmers[i][k] = '\0';
    }
    return ret;
};

/** Returns an mkmh_minimizer struct, equivalent to a tuple(kmer, position, kmer length), for every position in the genome **/
vector<mkmh_minimizer> Minimizers::kmer_tuples(string seq, int k)
{
    vector<string> kmers = this->kmerize(seq, k);
    vector<mkmh_minimizer> tups(kmers.size());
    for (int i = 0; i < kmers.size(); i++)
    {
        mkmh_minimizer mm;
        mm.seq = kmers[i];
        mm.pos = i;
        mm.length = k;
        tups[i] = mm;
    }

    return tups;
}

template <typename T>
vector<T> Minimizers::v_set(vector<T> kmers)
{
    set<T> s = set<T>(kmers.begin(), kmers.end());
    vector<T> ret = vector<T>(s.begin(), s.end());
    return ret;
}

/** Finds the (w, k) minimizers of a string **/
vector<mkmh_minimizer> Minimizers::getMinimizers(string &seq)
{
    vector<mkmh_minimizer> ret;
    vector<mkmh_minimizer> kmert = kmer_tuples(seq, this->k);
    int i = 0;
    for (i = 0; i + this->w < kmert.size(); ++i)
    {
        // get and sort kmers in window (i, i + w)
        vector<mkmh_minimizer> window_kmers(kmert.begin() + i, kmert.begin() + i + this->w);
        std::sort(window_kmers.begin(), window_kmers.end());
        ret.push_back(*(window_kmers.begin()));
    }
    return v_set(ret);
}

void Minimizers::extractKmers()
{
    vector<mkmh_minimizer> ret;
    vector<mkmh_minimizer> kmert;



    std::string id;
    std::string seq;

    for(unsigned long seq_num =0; seq_num < seqan::length(this->ids); seq_num++) {

        id = std::string((char *) seqan::toCString(this->ids[seq_num]));
        seq = std::string((char *) seqan::toCString(this->seqs[seq_num]));
        kmert = kmer_tuples(seq, this->k);

        for (unsigned long i = 0; i + this->w < kmert.size(); ++i)
        {
            // get and sort kmers in window (i, i + w)
            vector<mkmh_minimizer> window_kmers(kmert.begin() + i, kmert.begin() + i + this->w);
            std::sort(window_kmers.begin(), window_kmers.end());
            ret.push_back(*(window_kmers.begin()));
        }

        for (auto z : v_set(ret))
        {
            this->kmers[id].push_back(z.seq);
        }

    }

}