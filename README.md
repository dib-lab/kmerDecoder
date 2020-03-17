# kmerDecoder

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/675f273c446f45bebb5b8e354e24bccb)](https://app.codacy.com/app/mr-eyes/kmerDecoder?utm_source=github.com&utm_medium=referral&utm_content=mr-eyes/kmerDecoder&utm_campaign=Badge_Grade_Dashboard)
[![Build Status](https://travis-ci.org/dib-lab/kmerDecoder.svg?branch=master)](https://travis-ci.org/dib-lab/kmerDecoder)

## Quick Setup (using CMake)

Create `CMakeLists.txt` file in your project directory

```cmake
cmake_minimum_required (VERSION 3.4)
project (KD_Test C CXX)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fopenmp")

add_subdirectory(kmerDecoder)

add_executable (kdtest ${CPP_SOURCE_FILE})

target_link_libraries(kdtest kmerDecoder)

```

Use the generated static lib either using Cmake `add_subdirectory()` or by linking to compilation command.

## Usage Example

create `sample.fa`

```fasta
>SAMPLE
ACGTAGCATGCATGACGATGCTAGCGT
```

```cpp
#include "kmerDecoder.hpp"

void extract(kmerDecoder *KD);

int main() {
  std::string filename = "sample.fa";
  unsigned int chunk_size = 1;

  /*
  Create three kmerDecoder objects
  */

  // 1- Kmers Mode > Kmers(seqan_object, chunk_size, kSize)
  kmerDecoder *KMERS = new Kmers(filename, chunk_size, 15);

  // 2- Skipmers Mode > Skipmers(seqan_object, chunk_size, m, n, k)
  kmerDecoder *SKIPMERS = new Skipmers(filename, chunk_size, 2, 3, 10);

  // 3- Minimizers Mode > Minimizers(seqan_object, chunk_size, k, w)
  kmerDecoder *MINIMIZERS = new Minimizers(filename, chunk_size, 5, 10);

  // Start Extraction

  std::cout << "Mode: Kmers" << "\n";
  extract(KMERS);
  std::cout << "------------------------------------" << std::endl;

  std::cout << "Mode: Skipmers" << "\n";
  extract(SKIPMERS);
  std::cout << "------------------------------------" << std::endl;

  std::cout << "Mode: Minimizers" << "\n";
  extract(MINIMIZERS);
  std::cout << "------------------------------------" << std::endl;
}

void extract(kmerDecoder *KD) {
  while (!KD->end()) {
    KD->next_chunk();

    for (const auto &seq : *KD->getKmers()) {
      std::cout << "Read ID: " << seq.first << std::endl;
      for (const auto &kmer : seq.second) {
        std::cout << kmer.str << ": " << kmer.hash << std::endl;
      }
    }
  }
}
```

### Output

```text
Mode: Kmers
Read ID: sample
ACTGATCGTAGCTAG: 270607162
CTGATCGTAGCTAGC: 944645122
TGATCGTAGCTAGCA: 557632235
GATCGTAGCTAGCAT: 1027309074
ATCGTAGCTAGCATG: 937389716
TCGTAGCTAGCATGC: 386653658
CGTAGCTAGCATGCA: 472543720
GTAGCTAGCATGCAT: 769426497
TAGCTAGCATGCATC: 926137824
AGCTAGCATGCATCG: 58327299
GCTAGCATGCATCGT: 33644462
CTAGCATGCATCGTA: 873893451
TAGCATGCATCGTAG: 920819452
AGCATGCATCGTAGC: 343354516
GCATGCATCGTAGCT: 330181179
CATGCATCGTAGCTA: 94646772
ATGCATCGTAGCTAG: 1034175364
TGCATCGTAGCTAGC: 777708297
------------------------------------
Mode: Skipmers
Read ID: sample
ACGACGAGTA: 539404
CGACGAGTAC: 74282
GACGAGTACA: 374149
ACGAGTACAG: 210706
CGAGTACAGC: 854584
GAGTACAGCT: 1021218
AGTACAGCTC: 68631
GTACAGCTCT: 1045782
TACAGCTCTA: 75413
ACAGCTCTAC: 316100
CAGCTCTACT: 194024
AGCTCTACTG: 873565
GCTCTACTGC: 272874
CTATGTGCAG: 527346
TATGTGCAGA: 26575
ATGTGCAGAT: 191221
TGTGCAGATC: 455797
GTGCAGATCA: 334737
TGCAGATCAC: 32520
GCAGATCACG: 932927
CAGATCACGA: 532836
AGATCACGAG: 36901
GATCACGAGT: 471900
ATCACGAGTA: 638086
TCACGAGTAC: 467666
TGTCTACTGC: 1011830
GTCTACTGCT: 338374
TCTACTGCTG: 344123
CTACTGCTGA: 904801
TACTGCTGAT: 478035
ACTGCTGATG: 872943
CTGCTGATGT: 158651
TGCTGATGTG: 301791
GCTGATGTGC: 596234
CTGATGTGCA: 378027
TGATGTGCAG: 1019609
------------------------------------
Mode: Minimizers
Read ID: sample
ACTGA: 647
AGCAT: 485
AGCTA: 715
ATCGT: 25
```
