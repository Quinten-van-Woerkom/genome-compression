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

  auto file_size = std::filesystem::file_size(path);

  std::cout
    << "\n============================================================\n"
    << " Input\n"
    << "============================================================\n"
    << " Filename:                  " << path << '\n'
    << " Size:                      " << bytes_to_string(file_size) << '\n'
    << " Nucleotides (upper bound): " << file_size << '\n';

  auto start = std::chrono::high_resolution_clock::now();
  auto compressed = shared_tree{path};
  auto end = std::chrono::high_resolution_clock::now();
  auto construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  
  start = std::chrono::high_resolution_clock::now();
  compressed.sort_tree();
  end = std::chrono::high_resolution_clock::now();
  auto sorting_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  auto compressed_size = compressed.bytes();
  auto compressed_width = compressed.width();

  if (!histogram.empty()) {
    compressed.store_histogram(histogram);
  }

  auto new_path = path;
  new_path.replace_extension(".dag");
  compressed.save(new_path);

  std::cout
    << "\n============================================================\n"
    << " Output\n"
    << "============================================================\n"
    << " Filename:                  " << new_path << '\n'
    << " Size:                      " << bytes_to_string(compressed_size) << '\n'
    << " Nucleotides:               " << compressed_width*dna::size() << '\n'
    << " Compression ratio:         " << double(file_size)/double(compressed_size) << '\n';

  if (!histogram.empty())
    std::cout << " Histogram:                 " << histogram << '\n';

  std::cout
    << "\n============================================================\n"
    << " Tree dimensions\n"
    << "============================================================\n"
    << " Leaf size:                 " << dna::size() << " base pairs\n"
    << " Width:                     " << compressed_width << '\n'
    << " Depth:                     " << compressed.depth() << '\n'
    << " Leaves:                    " << compressed.leaf_count() << '\n'
    << " Nodes:                     " << compressed.node_count() << '\n';

  std::cout
    << "\n============================================================\n"
    << " Timings\n"
    << "============================================================\n"
    << " Tree construction:         " << construction_time.count() << " ms\n"
    << " Frequency sorting:         " << sorting_time.count() << " ms\n\n";

  // Short output for data collection purposes
  // std::cout << dna::size() << ',' << double(file_size)/double(compressed_size) << '\n';

  return 0;
}