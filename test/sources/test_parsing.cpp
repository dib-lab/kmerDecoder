#include <doctest/doctest.h>
#include <fstream>
#include "kmerDecoder.hpp"

// Test reading from file
string headers_filename = "test_data/mixed_headers.fa";

vector<string> fasta_to_headers_vec(string& filename) {
    std::vector<string> original_headers;
    string line;
    ifstream myfile(filename);

    if (myfile.is_open()) {
        while (getline(myfile, line)) {
            if (line.substr(0, 1) == ">") {
                original_headers.emplace_back(line.substr(1));
            }
        }
        myfile.close();
    }
    else {
        throw "Unable to open file";
    }

    return original_headers;
}

int number_of_seqs = 0;

int chunk_size = 1;
int kSize = 11;
kmerDecoder* PROT_KMERS = new aaKmers(headers_filename, chunk_size, kSize);

auto original_headers = fasta_to_headers_vec(headers_filename);
int golden_seq_count = original_headers.size();
int parsed_seqs_count = 0;
int seq_idx = -1;

TEST_CASE("Sequences Parsing test") {
    SUBCASE("Parsing headers") {
        while (!PROT_KMERS->end()) {
            PROT_KMERS->next_chunk();

            for (const auto& seq : *PROT_KMERS->getKmers()) {
                parsed_seqs_count++;
                CHECK_EQ(seq.first, original_headers[++seq_idx]);
            }
        }
    }

    SUBCASE("Parsed sequences number") {
        CHECK_EQ(parsed_seqs_count, golden_seq_count);
    }


}
