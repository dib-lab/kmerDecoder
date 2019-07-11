# kmerDecoder

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/675f273c446f45bebb5b8e354e24bccb)](https://app.codacy.com/app/mr-eyes/kmerDecoder?utm_source=github.com&utm_medium=referral&utm_content=mr-eyes/kmerDecoder&utm_campaign=Badge_Grade_Dashboard)
[![Build Status](https://travis-ci.org/mr-eyes/kmerDecoder.svg?branch=master)](https://travis-ci.org/mr-eyes/kmerDecoder)

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

  // Create three Seqan objects for the same file
  seqan::SeqFileIn SeqIn_kmers;
  seqan::SeqFileIn SeqIn_skipmers;
  seqan::SeqFileIn SeqIn_minimizers;

  // Initialize Seqan for the three modes
  seqan::open(SeqIn_kmers, seqan::toCString(filename));
  seqan::open(SeqIn_skipmers, seqan::toCString(filename));
  seqan::open(SeqIn_minimizers, seqan::toCString(filename));

  /*
  Create three kmerDecoder objects
  */

  // 1- Kmers Mode > Kmers(seqan_object, chunk_size, kSize)
  kmerDecoder *KMERS = new Kmers(SeqIn_kmers, chunk_size, 15);

  // 2- Skipmers Mode > Skipmers(seqan_object, chunk_size, m, n, k)
  kmerDecoder *SKIPMERS = new Skipmers(SeqIn_skipmers, chunk_size, 2, 3, 10);

  // 3- Minimizers Mode > Minimizers(seqan_object, chunk_size, k, w)
  kmerDecoder *MINIMIZERS = new Minimizers(SeqIn_minimizers, chunk_size, 5, 10);

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
        std::cout << kmer << std::endl;
      }
    }
  }
}
```
