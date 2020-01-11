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

  auto data = dna::random();
  auto leaf = pointer{data};
  if (!leaf.is_leaf()) {
    std::cerr << "<Tree node> Test failed: pointer to data should be a leaf node\n";
    ++errors;
  }

  if (leaf.data() != data) {
    std::cerr << "<Tree node> Test failed: pointer{data}.leaf() != data: " << leaf.data() << "!=" << data << '\n';
    ++errors;
  }

  auto node_pointer = pointer{1};
  if (node_pointer.is_leaf()) {
    std::cerr << "<Tree node> Test failed: pointer to points should not be a leaf node\n";
    ++errors;
  }

  auto empty = pointer{nullptr};
  if (!empty.empty()) {
    std::cerr << "<Tree node> Test failed: nullptr-initialized pointer should be empty\n";
    ++errors;
  }

  if (!errors) std::cout << "<Tree node> Finished without errors\n";

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
  auto buffered = fasta_reader{path};
  auto direct = std::vector<std::array<char, dna::size()>>(size/dna::size());
  auto file = std::ifstream{path, std::ios::binary};
  auto errors = 0;

  file.read(reinterpret_cast<char*>(direct.data()), size);

  if (direct.size() != size/dna::size()) {
    std::cerr << "<File reader> Test failed: file path size does not match read size within accuracy bounds\n"
      << "File path: " << size << '\n'
      << "Direct read: " << direct.size()*dna::size() << '\n';
    ++errors;
  }

  auto i = 0u;
  for (auto b : buffered) {
    if (b != std::string_view{direct[i].data(), direct[i].size()}) {
      std::cerr << "<File reader> Test failed: buffered[i] != direct[i] for i = " << i << " out of " << direct.size() - 1 << '\n'
        << "\tbuffered[i]\t= " << b << '\n'
        << "\tdirect[i]\t= " << std::string_view{direct[i].data(), direct[i].size()} << '\n';
      ++errors;
    }
    ++i;
  }

  if (i != size/dna::size()) {
    std::cerr << "<File reader> Test failed: file path size does not match buffered read size within accuracy bounds\n"
      << "File path: " << size << '\n'
      << "Buffered read: " << i*dna::size() << '\n';
    ++errors;
  }

  if (!errors) std::cout << "<File reader> Finished without errors\n";

  return errors;
}

auto test_tree_factory() -> int {
  auto path = "data/chmpxx";
  auto data = read_genome(path);
  auto compressed = shared_tree::create_balanced(data);
  auto errors = 0;

  if (data.size() != compressed.width()) {
    std::cerr << "<Tree factory> Test failed: raw data size (" << data.size()
      << ") does not match compressed data size (" << compressed.width() << ", " << compressed.node_count() << " nodes)\n";
    ++errors;
  }

  auto i = 0;
  for (const auto c : compressed) {
    if (data[i] != c) {
      std::cerr << "<Tree factory> Test failed: data[i] != compressed[i] for i = " << i << " out of " << data.size() - 1 << '\n'
        << "\tdata[i]\t\t= " << data[i] << '\n'
        << "\tcompressed[i]\t= " << c << '\n';
      ++errors;
    }
    ++i;
  }

  if (!errors) std::cout << "<Tree factory> Finished without errors\n";

  return errors;
}

int main(int argc, char* argv[]) {
  auto errors = test_chunks() + test_file_reader() + test_node() + test_tree_factory();
  if (errors) std::cerr << "Not all tests passed\n";
  return errors;
}