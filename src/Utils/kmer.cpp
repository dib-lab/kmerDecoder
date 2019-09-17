
#include <stdio.h>
#include <string>
#include "Utils/kmer.h"

using namespace std;


/**
* Converts a string of "ATCG" to a uint64_t
* where each character is represented by using only two bits
*/
uint64_t kmer::str_to_int(string str)
{
	uint64_t strint = 0;
	for (auto it = str.begin(); it != str.end(); it++) {
		uint8_t curr = 0;
		switch (*it) {
			case 'A': { curr = DNA_MAP::A; break; }
			case 'T': { curr = DNA_MAP::T; break; }
			case 'C': { curr = DNA_MAP::C; break; }
			case 'G': { curr = DNA_MAP::G; break; }
		}
		strint = strint | curr;
		strint = strint << 2;
	}
	return strint >> 2;
}
uint64_t kmer::str_to_canonical_int(string str)
{
	uint64_t kmerI=kmer::str_to_int(str);
	uint64_t kmerIR=kmer::reverse_complement(kmerI,str.size());
	if (kmer::compare_kmers(kmerI, kmerIR))
		return kmerI;
	else
		return kmerIR;
}
 string kmer::canonicalKmer(string kmer){
	uint64_t kmerI = kmer::str_to_int(kmer);
	uint64_t kmerIR = kmer::reverse_complement(kmerI, kmer.size());
	uint64_t item;
	if (kmer::compare_kmers(kmerI, kmerIR))
		return kmer;
	else
		return kmer::int_to_str(kmerIR,kmer.size());
}
/**
* Converts a uint64_t to a string of "ACTG"
* where each character is represented by using only two bits
*/
string kmer::int_to_str(uint64_t kmer, uint32_t K)
{
	uint8_t base;
	string str;
	for (int i=K; i>0; i--) {
		base = (kmer >> (i*2-2)) & 3ULL;
		char chr = kmer::map_int(base);
		str.push_back(chr);
	}
	return str;
}

/* Calculate the revsese complement of a kmer */
uint64_t kmer::reverse_complement(uint64_t kmer, uint32_t K)
{
	uint64_t rc = 0;
	uint8_t base = 0;
	for (int i=0; i<K; i++) {
		base = kmer & 3ULL;
		base = reverse_complement_base(base);
		kmer >>= 2;
		rc |= base;
		rc <<= 2;
	}
	rc >>=2;
	return rc;
}
