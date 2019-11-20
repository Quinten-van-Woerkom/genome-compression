/**
 *  Unit tests for the implementation.
 */

#include <cassert>
#include <iostream>

#include "dna.h"
#include "shared_tree.h"

bool test_tree() {
  auto path = "data/chmpxx";
  auto data = read_genome(path);
  auto compressed = shared_tree<dna>{data};

  // This could be much more efficient.
  for (auto i = 0ul; i < data.size(); ++i) {
    if (data[i] != compressed[i]) return 1;
  }
  return 0;
}

int main(int argc, char* argv[]) {
  auto errors = test_tree();
  if (errors) std::cerr << "Not all tests passed\n";
  return errors;
}