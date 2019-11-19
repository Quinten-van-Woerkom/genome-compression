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
  // Todo: provide buffered file reader implementation

  auto data = read_genome(path);
  std::cout << "Data length: " << data.size() << " elements\n";
  auto compressed = shared_binary_tree<dna>::construct_from(data);
  std::cout << "Done compressing\n"
    << "Final size: " << compressed.size() << " nodes\n"
    << "Data length: " << compressed.length() << " elements\n";
}