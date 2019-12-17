/**
 *  Unit tests for the implementation.
 */

#include <array>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "dna.h"
#include "file_reader.h"
#include "shared_tree.h"
#include "utility.h"

auto test_node() -> int {
  auto errors = 0;

  auto leaf = node::pointer{dna::random()};
  if (!leaf.is_leaf()) {
    std::cerr << "<Tree node> Test failed: pointer to data should be a leaf node\n";
    ++errors;
  }

  auto inner_node = node{leaf, leaf};
  auto node_pointer = node::pointer{inner_node};
  if (node_pointer.is_leaf()) {
    std::cerr << "<Tree node> Test failed: pointer to points should not be a leaf node\n";
    ++errors;
  }

  auto empty = node::pointer{nullptr};
  if (!empty.empty()) {
    std::cerr << "<Tree node> Test failed: nullptr-initialized pointer should be empty\n";
    ++errors;
  }

  if (!errors) std::cout << "<Tree node> Finished without errors\n";

  return errors;
}

auto test_file_reader() -> int {
  auto path = "data/hehcmv";
  auto size = std::filesystem::file_size(path);
  auto buffered = fasta_reader{path, 32768, dna::size()};
  auto direct = std::vector<std::array<char, dna::size()>>(size);
  auto file = std::ifstream{path, std::ios::binary};
  auto errors = 0;

  file.read(reinterpret_cast<char*>(direct.data()), direct.size());

  auto i = 0;
  for (const auto& [b, d] : zip(buffered, direct)) {
    if (b != std::string_view{d.data(), d.size()}) {
      std::cerr << "<File reader> Test failed: buffered[i] != direct[i] for i = " << i << " out of " << direct.size() - 1 << '\n'
        << "\tbuffered[i]\t= " << b << '\n'
        << "\tdirect[i]\t= " << std::string_view{d.data(), d.size()} << '\n';
      ++errors;
    }
    ++i;
  }

  if (!errors) std::cout << "<File reader> Finished without errors\n";

  return errors;
}

auto test_tree_layerwise() -> int {
  auto path = "data/hehcmv";
  auto data = read_genome(path);
  auto compressed = shared_tree::create_balanced_layerwise(data);
  auto errors = 0;

  if (data.size() != compressed.length()) {
    std::cerr << "<Tree layerwise> Test failed: raw data size (" << data.size()
      << ") does not match compressed data size (" << compressed.length() << ")\n";
    ++errors;
  }

  auto i = 0;
  for (const auto& [d, c] : zip(data, compressed)) {
    if (d != c) {
      std::cerr << "<Tree layerwise> Test failed: data[i] != compressed[i] for i = " << i << " out of " << data.size() - 1 << '\n'
        << "\tdata[i]\t\t= " << d << '\n'
        << "\tcompressed[i]\t= " << c << '\n';
      ++errors;
    }
    ++i;
  }

  if (!errors) std::cout << "<Tree layerwise> Finished without errors\n";

  return errors;
}

auto test_tree_pairwise() -> int {
  auto path = "data/hehcmv";
  auto data = read_genome(path);
  auto compressed = shared_tree::create_balanced_pairwise(data);
  auto errors = 0;

  if (data.size() != compressed.length()) {
    std::cerr << "<Tree pairwise> Test failed: raw data size (" << data.size()
      << ") does not match compressed data size (" << compressed.length() << ")\n";
    ++errors;
  }

  auto i = 0;
  for (const auto& [d, c] : zip(data, compressed)) {
    if (d != c) {
      std::cerr << "<Tree pairwise> Test failed: data[i] != compressed[i] for i = " << i << " out of " << data.size() - 1 << '\n'
        << "\tdata[i]\t\t= " << d << '\n'
        << "\tcompressed[i]\t= " << c << '\n';
      ++errors;
    }
    ++i;
  }

  if (!errors) std::cout << "<Tree pairwise> Finished without errors\n";

  return errors;
}

auto test_tree_factories() -> int {
  auto path = "data/hehcmv";
  auto data = read_genome(path);
  auto layerwise = shared_tree::create_balanced_layerwise(data);
  auto pairwise = shared_tree::create_balanced_pairwise(data);
  auto errors = 0;

  if (layerwise.length() != pairwise.length()) {
    std::cerr << "<Tree factories> Test failed: layerwise-constructed tree length (" << layerwise.length()
      << ") does not match pairwise-constructed tree length (" << pairwise.length() << ")\n";
    ++errors;
  }

  auto i = 0;
  for (const auto& [l, p] : zip(layerwise, pairwise)) {
    if (l != p) {
      std::cerr << "<Tree factories> Test failed: layerwise[i] != pairwise[i] for i = " << i << " out of " << layerwise.size() - 1 << '\n'
        << "\tlayerwise[i]\t= " << l << '\n'
        << "\tpairwise[i]\t= " << p << '\n';
      ++errors;
    }
    ++i;
  }

  if (!errors) std::cout << "<Tree factories> Finished without errors\n";

  return errors;
}

int main(int argc, char* argv[]) {
  auto errors = test_tree_factories() + test_tree_layerwise() + test_tree_pairwise() + test_file_reader() + test_node();
  if (errors) std::cerr << "Not all tests passed\n";
  return errors;
}