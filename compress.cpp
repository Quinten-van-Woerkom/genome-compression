/**
 *  Entry point to the application.
 *  Compresses the given file using a canonicalized directed acyclic graph.
 */

#include <filesystem>
#include <iostream>
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

  dna left{"ACGTTGCA"};
  dna right{"TGACTGAC"};
  node<dna> dna_node{left, right};
  std::cout << "Hello, world!\n";
  std::cout << dna_node.right_leaf() << '\n';
}