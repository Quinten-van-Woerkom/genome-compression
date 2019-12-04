/**
 *  Unit tests for the implementation.
 */

#include <array>
#include <cassert>
#include <iostream>

#include "dna.h"
#include "shared_tree.h"
#include "utility.h"

auto test_node() -> int {
  auto errors = 0;

  auto leaf = node::pointer{dna::random()};
  if (!leaf.is_leaf()) {
    std::cerr << "Test failed: pointer to data should be a leaf node\n";
    ++errors;
  }

  auto inner_node = node{leaf, leaf};
  auto node_pointer = node::pointer{inner_node};
  if (node_pointer.is_leaf()) {
    std::cerr << "Test failed: pointer to points should not be a leaf node\n";
    ++errors;
  }

  auto empty = node::pointer{nullptr};
  if (!empty.empty()) {
    std::cerr << "Test failed: nullptr-initialized pointer should be empty\n";
    ++errors;
  }

  return errors;
}

auto test_tree() -> int {
  auto path = "data/chmpxx";
  auto data = read_genome(path);
  auto compressed = shared_tree{data};
  auto errors = 0;

  if (data.size() != compressed.length()) {
    std::cerr << "Test failed: raw data size (" << data.size()
      << ") does not match compressed data size (" << compressed.length() << ")\n";
    ++errors;
  }

  auto i = 0;
  for (const auto& [d, c] : zip(data, compressed)) {
    ++i;
    if (d != c) {
      std::cerr << "Test failed: data[i] != compressed[i] for i = " << i << " out of " << data.size() - 1 << '\n'
        << "\tdata[i]\t\t= " << d << '\n'
        << "\tcompressed[i]\t= " << c << '\n';
      ++errors;
    }
  }

  return errors;
}

int main(int argc, char* argv[]) {
  auto errors = test_tree() + test_node();
  if (errors) std::cerr << "Not all tests passed\n";
  return errors;
}