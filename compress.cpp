/**
 *  Entry point to the application.
 *  Compresses the given file using a canonicalized directed acyclic graph.
 */

#include <filesystem>
#include <iostream>
#include "dna.h"

int main(int argc, char* argv[]) {
  // std::string path = "data/chmpxx";
  // if (argc > 1) {
  //     path = argv[1];
  // }
  // std::cout << "Compressing " << path << '\n';
  // Todo: provide buffered file reader implementation

  dna test{"ACGTTGCA"};
  std::cout << test << '\n';
}