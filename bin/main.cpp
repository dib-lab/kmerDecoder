#include "kmerDecoder.hpp"

void extract(kmerDecoder *KD);

int main() {
  std::string filename = "sample.fa";
  std::string prot_filename = "aaSample.fa";
  unsigned int chunk_size = 100;

  /*
  Create three kmerDecoder objects
  */

  // 1- Kmers Mode > Kmers(file_name, chunk_size, kSize)
  kmerDecoder *KMERS = new Kmers(filename, chunk_size, 15);

  // 2- Skipmers Mode > Skipmers(file_name, chunk_size, m, n, k)
  kmerDecoder *SKIPMERS = new Skipmers(filename, chunk_size, 2, 3, 10);

  // 3- Minimizers Mode > Minimizers(file_name, chunk_size, k, w)
  kmerDecoder *MINIMIZERS = new Minimizers(filename, chunk_size, 5, 10);

   // 4- Proteins Mode > aaKmers(file_name, chunk_size, k, w)
  kmerDecoder *PROT_KMERS = new aaKmers(prot_filename, chunk_size, 11);

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

  std::cout << "Mode: Protein Kmers" << "\n";
  extract(PROT_KMERS);
  std::cout << "------------------------------------" << std::endl;

  // Freeing the memory
    delete KMERS;
    delete SKIPMERS;
    delete MINIMIZERS;
    delete PROT_KMERS;


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