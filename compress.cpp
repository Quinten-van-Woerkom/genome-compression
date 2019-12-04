/**
 *  Entry point to the application.
 *  Compresses the given file using a canonicalized directed acyclic graph.
 */

#include <chrono>
#include <filesystem>
#include <iostream>
#include <vector>

#include "dna.h"
#include "shared_tree.h"
#include "file_reader.h"

int main(int argc, char* argv[]) {
  std::string path = "data/chmpxx";
  if (argc > 1) {
      path = argv[1];
  }
  std::cout << "Compressing " << path << "...";

  auto start = std::chrono::high_resolution_clock::now();

  auto data = read_genome(path);
  // auto data = fasta_reader{path};
  auto compressed = shared_tree{data};

  auto end = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  auto compressed_size = compressed.size()*sizeof(node);
  auto file_size = std::filesystem::file_size(path);

  std::cout << " Done!\n"
    << "Number of base pairs: " << compressed.length()*dna::size() << '\n'
    << "File size: " << file_size << " bytes\n"
    << "Compressed size: " << compressed.size() << " nodes (" << compressed_size << " bytes)\n"
    << "Compression ratio: " << double(compressed_size)/double(file_size) << '\n'
    << "Time: " << time.count() << "ms\n";

  return 0;
}