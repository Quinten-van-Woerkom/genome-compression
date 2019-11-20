/**
 *  Entry point to the application.
 *  Compresses the given file using a canonicalized directed acyclic graph.
 */

#include <filesystem>
#include <iostream>
#include <vector>

#include "dna.h"
#include "shared_tree.h"

int main(int argc, char* argv[]) {
  std::string path = "data/chmpxx";
  if (argc > 1) {
      path = argv[1];
  }
  std::cout << "Compressing " << path << '\n';

  auto data = read_genome(path);
  auto compressed = shared_binary_tree<dna>{data};

  auto compressed_size = compressed.size()*sizeof(decltype(compressed)::node_type);
  auto uncompressed_size = compressed.length()*sizeof(decltype(compressed)::value_type);

  std::cout << "Done compressing\n"
    << "Final size: " << compressed.size() << " nodes (" << compressed_size << " bytes)\n"
    << "Data length: " << compressed.length() << " elements (" << uncompressed_size << " bytes)\n"
    << "Compression ratio: " << double(compressed_size)/double(uncompressed_size) << '\n';
}