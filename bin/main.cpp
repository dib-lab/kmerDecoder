#include "kmerDecoder.hpp"

void extract(kmerDecoder *KD);

int main() {
  std::string filename = "sample.fa";
  unsigned int chunk_size = 100;

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

  // Freeing the memory
    delete KMERS;
    delete SKIPMERS;
    delete MINIMIZERS;

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