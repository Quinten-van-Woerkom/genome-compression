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
  auto compressed = shared_tree{data};

  // This could be much more efficient.
  auto errors = 0;
  for (auto i = 0ul; i < data.size(); ++i) {
    if (data[i] != compressed[i]) {
      std::cerr << "Test failed: data[i] != compressed[i] for i = " << i << " out of " << data.size() - 1 << '\n'
        << "\tdata[i]\t\t= " << data[i] << '\n'
        << "\tcompressed[i]\t= " << compressed[i] << '\n';
      ++errors;
    }
  }
  return errors;
}

int main(int argc, char* argv[]) {
  auto errors = test_tree();
  if (errors) std::cerr << "Not all tests passed\n";
  return errors;
}