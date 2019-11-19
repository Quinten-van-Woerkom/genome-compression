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
  // std::string path = "data/chmpxx";
  // if (argc > 1) {
  //     path = argv[1];
  // }
  // std::cout << "Compressing " << path << '\n';
  // Todo: provide buffered file reader implementation

  dna test{"ACGTTGCA"};
  std::cout << test << '\n';

  auto data = std::vector<dna>{
    {"ACGTTGCA"}, {"TGACTGAC"}, {"TTTTAAAA"}
  };

  auto tree = shared_binary_tree<dna>::construct_from(data);
}