/**
 *  Entry point to the application.
 *  Compresses the given file using a canonicalized directed acyclic graph.
 */

#include <chrono>
#include <filesystem>
#include <iostream>
#include <vector>

#include "balanced_shared_tree.h"
#include "dna.h"
#include "fasta_reader.h"

int main(int argc, char* argv[]) {
  std::string path = "data/chmpxx";
  std::string histogram;
  if (argc > 1) path = argv[1];
  if (argc > 2) histogram = argv[2];

  std::cout << "Compressing " << path << "...\n";
  std::cout.flush();

  auto start = std::chrono::high_resolution_clock::now();

  // auto data = fasta_reader{path};
  // auto compressed = shared_tree::create_balanced(data);
  auto compressed = balanced_shared_tree{path};

  auto end = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  std::cout << "Done!\n" << std::flush;

  auto compressed_size = compressed.node_count();
  auto compressed_width = compressed.width();
  auto file_size = std::filesystem::file_size(path);
  if (!histogram.empty()) {
    compressed.histogram(histogram);
    std::cout << "Stored node reference histogram at " << histogram << '\n';
  }

  std::cout
    << "Number of data units: " << compressed_width << '\n'
    << "Number of base pairs: " << compressed_width*dna::size() << '\n'
    << "File size: " << file_size << " bytes\n"
    << "Compressed size: " << compressed_size*sizeof(balanced_shared_tree::node) << " bytes (" << compressed_size << " nodes)\n"
    << "Compression ratio: " << double(compressed_size*sizeof(balanced_shared_tree::node))/double(file_size) << '\n'
    << "Time: " << time.count() << "ms\n";

  return 0;
}