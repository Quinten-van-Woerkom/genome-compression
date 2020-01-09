/**
 *  Unit tests for the implementation.
 */

#include <array>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "dna.h"
#include "fasta_reader.h"
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

auto test_zip() -> int {
  auto a = std::array{1, 2, 3};
  auto b = std::array{4, 5, 6};
  auto errors = 0;
  auto index = 1;

  for (auto [ae, be] : zip(a, b)) {
    if (ae != index) {
      std::cerr << "Element mismatch in a: " << ae << " != " << index << '\n';
      ++errors;
    }
    if (be != index+3) {
      std::cerr << "Element mismatch in b: " << be << " != " << index+3 << '\n';
      ++errors;
    }
    ++index;

    ae = 0;
  }

  for (const auto& ae : a) {
    if (ae != 0) {
      std::cerr << "Element mismatch in a: " << ae << " != " << 0 << '\n';
      ++errors;
    }
  }
  
  if (!errors) std::cout << "<Zip> Finished without errors\n";

  return errors;
}

auto test_chunks() -> int {
  auto a = std::array{1, 2, 3, 4, 5, 6, 7, 8};
  auto errors = 0;
  auto index = 1;

  for (auto chunk : chunks(a, 2)) {
    for (auto v : chunk) {
      if (v != index) {
        std::cerr << "<Chunks> Value mismatch: " << v << " != " << index << '\n';
        ++errors;
      }
      ++index;
    }
  }
  
  if (!errors) std::cout << "<Chunks> Finished without errors\n";

  return errors;
}

auto test_file_reader() -> int {
  auto path = "data/chmpxx";
  auto size = std::filesystem::file_size(path);
  auto buffered = fasta_reader{path, dna::size()};
  auto direct = std::vector<std::array<char, dna::size()>>(size/dna::size());
  auto file = std::ifstream{path, std::ios::binary};
  auto errors = 0;

  file.read(reinterpret_cast<char*>(direct.data()), size);

  auto i = 0;
  for (auto b : buffered) {
    if (b != std::string_view{direct[i].data(), direct[i].size()}) {
      std::cerr << "<File reader> Test failed: buffered[i] != direct[i] for i = " << i << " out of " << direct.size() - 1 << '\n'
        << "\tbuffered[i]\t= " << b << '\n'
        << "\tdirect[i]\t= " << std::string_view{direct[i].data(), direct[i].size()} << '\n';
      ++errors;
    }
    ++i;
  }

  if (!errors) std::cout << "<File reader> Finished without errors\n";

  return errors;
}

auto test_tree() -> int {
  auto left = dna::random();
  auto right = dna::random();
  auto single = dna::random();
  auto left_node = node{left, right};
  auto right_node = node{single};
  auto root = node{left_node, right_node};
  auto errors = 0;

  if (root[0] != left) {
    std::cerr << "<Tree> Test failed: root[0] != left\n"
      << "\troot[0] =\t" << root[0] << '\n'
      << "\tleft =\t" << left << '\n';
    ++errors;
  }

  if (root[1] != right) {
    std::cerr << "<Tree> Test failed: root[1] != right\n"
      << "\troot[1] =\t" << root[1] << '\n'
      << "\tright =\t" << right << '\n';
    ++errors;
  }

  if (root[2] != single) {
    std::cerr << "<Tree> Test failed: root[2] != single\n"
      << "\troot[2] =\t" << root[2] << '\n'
      << "\tsingle =\t" << single << '\n';
    ++errors;
  }

  if (!errors) std::cout << "<Tree> Finished without errors\n";
  return errors;
}

auto test_tree_factory() -> int {
  auto path = "data/chmpxx";
  auto data = read_genome(path);
  auto compressed = shared_tree::create_balanced(data);
  auto errors = 0;

  if (data.size() != compressed.width()) {
    std::cerr << "<Tree factory> Test failed: raw data size (" << data.size()
      << ") does not match compressed data size (" << compressed.width() << ")\n";
    ++errors;
  }

  auto i = 0;
  for (const auto [d, c] : zip(data, compressed)) {
    if (d != c) {
      std::cerr << "<Tree factory> Test failed: data[i] != compressed[i] for i = " << i << " out of " << data.size() - 1 << '\n'
        << "\tdata[i]\t\t= " << d << '\n'
        << "\tcompressed[i]\t= " << c << '\n';
      ++errors;
    }
    ++i;
  }

  if (!errors) std::cout << "<Tree factory> Finished without errors\n!";

  return errors;
}

int main(int argc, char* argv[]) {
  auto errors = test_zip() + test_chunks() + test_file_reader() + test_node() + test_tree() + test_tree_factory();
  if (errors) std::cerr << "Not all tests passed\n";
  return errors;
}