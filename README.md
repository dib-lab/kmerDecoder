# kmerDecoder

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/b77991fd08eb4a98a05ac9ea0b812753)](https://app.codacy.com/app/mr-eyes/kmerDecoder?utm_source=github.com&utm_medium=referral&utm_content=mr-eyes/kmerDecoder&utm_campaign=Badge_Grade_Settings)
[![Build Status](https://travis-ci.org/mr-eyes/kmerDecoder.svg?branch=master)](https://travis-ci.org/mr-eyes/kmerDecoder)


## Build

Use the generated static lib either using Cmake `add_subdirectory()` or by linking to compilation command.




## Usage Example

```cpp

#include "inputModule.hpp"

std::string seq = "ACGTAGCATGCATGACGATGCTAGCGT";

// Kmers Parameters
int kSize = 15; // kmer size

// Skipmers Parameters
int sk_m = 2;
int sk_n = 3;
int sk_k = 10

// Minimizers Paramters
int min_k = 5
int min_w = 10

```

### Extract Kmers

```cpp

KDecoder *KMERS = new Kmers(kSize);

for(const auto &kmer : *KMERS->getKmers(seq))
    std::cout << kmer << std::endl;
```
Output:
```shell=
ACGTAGCATGCATGA
CGTAGCATGCATGAC
GTAGCATGCATGACG
TAGCATGCATGACGA
AGCATGCATGACGAT
GCATGCATGACGATG
CATGCATGACGATGC
ATGCATGACGATGCT
TGCATGACGATGCTA
GCATGACGATGCTAG
CATGACGATGCTAGC
ATGACGATGCTAGCG
TGACGATGCTAGCGT
```

### Extract Skipmers

```cpp
KDecoder *SKIPMERS = new Skipmers(sk_m,sk_n,sk_k);

for(const auto &kmer : *SKIPMERS->getKmers(seq))
    std::cout << kmer << std::endl;
```
Output:
```shell=
ACTACAGCTG
CTACAGCTGC
TACAGCTGCG
ACAGCTGCGT
CAGCTGCGTG
AGCTGCGTGT
GCTGCGTGTA
CTGCGTGTAC
TGCGTGTACG
CGAGATCAGA
GAGATCAGAG
AGATCAGAGA
GATCAGAGAG
ATCAGAGAGC
TCAGAGAGCA
CAGAGAGCAG
AGAGAGCAGG
GAGAGCAGGT
GTGCTGATAC
TGCTGATACA
GCTGATACAT
CTGATACATC
TGATACATCT
GATACATCTG
ATACATCTGC
TACATCTGCT
```

### Extract Minimizers

```cpp
KDecoder *MINIMZERS = new Minimzers(min_k,min_w);

for(const auto &kmer : *MINIMZERS->getKmers(seq))
    std::cout << kmer << std::endl;
```
Output:
```
ACGAT
ACGTA
AGCAT
```
