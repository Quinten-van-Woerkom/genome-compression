/**
 *  Entry point to the application.
 *  Compresses the given file using a canonicalized directed acyclic graph.
 */

#include <chrono>
#include <filesystem>
#include <iostream>
#include <vector>

#include "shared_tree.h"
#include "dna.h"
#include "fasta_reader.h"

int main(int argc, char* argv[]) {
  std::filesystem::path path = "data/chmpxx";
  std::filesystem::path histogram;
  if (argc > 1) path = argv[1];
  if (argc > 2) histogram = argv[2];

  std::cout << "Compressing " << path << "...\n";
  std::cout.flush();

  auto start = std::chrono::high_resolution_clock::now();

  // auto data = fasta_reader{path};
  // auto compressed = shared_tree::create_balanced(data);
  auto compressed = shared_tree{path};

  auto end = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  std::cout << "Done!\n" << std::flush;

  std::cout << "Sorting pointers based on frequency\n";
  compressed.sort_tree();

  auto compressed_size = compressed.bytes();
  auto compressed_width = compressed.width();
  auto file_size = std::filesystem::file_size(path);
  if (!histogram.empty()) {
    compressed.store_histogram(histogram);
    std::cout << "Stored node reference histogram at " << histogram << '\n';
  }

  path.replace_extension(".dag");
  compressed.save(path);
  std::cout << "Stored compressed sequence at " << path << '\n';

  std::cout
    << "Number of data units: " << compressed_width << '\n'
    << "Number of base pairs: " << compressed_width*dna::size() << '\n'
    << "File size: " << file_size << " bytes\n"
    << "Compressed size: " << compressed_size << " bytes (" << compressed.node_count() << " nodes, " << compressed.leaf_count() << " leaves)\n"
    << "Compression ratio: " << double(compressed_size)/double(file_size) << '\n'
    << "Time: " << time.count() << "ms\n"
    << "Last element: " << compressed[compressed_width-1] << '\n';

  return 0;
}