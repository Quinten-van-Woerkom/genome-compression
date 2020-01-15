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

auto test_pointer() -> int {
  auto errors = 0;
  auto leaf = dna{"ACGTACGTACGTACG"};
  auto basis = pointer{leaf, false, false};
  auto transposed = basis.transposed();
  auto mirrored = basis.mirrored();

  if (basis == transposed) {
    std::cerr << "<Tree pointer> Pointers can never match when transposed\n";
    ++errors;
  }

  if (basis == mirrored) {
    std::cerr << "<Tree pointer> Pointers that are not invariant under mirroring should not match when mirrored\n";
    ++errors;
  }

  if (!errors) std::cout << "<Tree pointer> Finished without errors\n";
  return errors;
}

auto test_node() -> int {
  auto errors = 0;

  auto data = dna::random();
  auto leaf = pointer{data, false, false};
  if (!leaf.is_leaf()) {
    std::cerr << "<Tree node> Test failed: pointer to data should be a leaf node\n";
    ++errors;
  }

  if (leaf.leaf() != data) {
    std::cerr << "<Tree node> Test failed: pointer{data}.leaf() != data: " << leaf.leaf() << "!=" << data << '\n';
    ++errors;
  }

  auto node_pointer = pointer{1, false, false};
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

auto test_tree_similarity_transforms() -> int {
  auto errors = 0;
  auto data = dna{"AAAAAAAAAGAAAAC"};
  auto other = dna{"AAAAAAAAAAAAAAA"};
  auto basis = data;
  auto transposed = data.transposed();
  auto mirrored = data.mirrored();
  auto both = data.transposed().mirrored();
  auto layer = std::vector<pointer>{};

  auto tree = shared_tree{};
  tree.emplace_node(layer, basis, transposed);
  tree.emplace_node(layer, both, mirrored);
  tree.emplace_node(layer, layer[0], layer[1]);

  auto basis_node = node{basis.canonical(), basis.canonical()};
  auto mirrored_node = node{mirrored.canonical(), mirrored.canonical()};
  auto transposed_node = node{transposed.canonical(), transposed.canonical()};
  auto mix_node = node{basis.canonical(), transposed.canonical()};
  auto mirrored_mix_node = node{both.canonical(), mirrored.canonical()};

  if (basis_node != mirrored_node
    && std::hash<node>()(basis_node) == std::hash<node>()(mirrored_node)) {
    std::cerr << "<Tree similarity> Mirrored nodes should compare to be equal\n";
    ++errors;
  }

  if (basis_node != transposed_node
    && std::hash<node>()(basis_node) == std::hash<node>()(transposed_node)) {
    std::cerr << "<Tree similarity> Transposed nodes should compare to be equal\n";
    ++errors;
  }

  if (mix_node != mirrored_mix_node
    && std::hash<node>()(mix_node) == std::hash<node>()(mirrored_mix_node)) {
    std::cerr << "<Tree similarity> Transposed and mirrored nodes should compare to be equal\n";
    ++errors;
  }

  auto left = node{basis.canonical(), other.canonical()};
  auto right = node{other.canonical(), basis.canonical()};
  if (left == right || std::hash<node>()(left) == std::hash<node>()(right)) {
    std::cerr << "<Tree similarity> Nodes with reversed children that are not mirrored should not match if they are not invariant under transformation\n";
    ++errors;
  }

  if (tree.node_count() != 2) {
    std::cerr << "<Tree similarity> Similar nodes do not merge: node count is "
      << tree.node_count() << ", exceeding the expected 2\n";
    ++errors;
  }

  if (!errors) std::cout << "<Tree similarity> Finished without errors\n";
  return errors;
}

auto test_tree_transposition() -> int {
  auto errors = 0;
  auto data = std::vector<dna>{};
  data.emplace_back("AAAAAAAAAAAAAAA");
  data.emplace_back("AAAAAAAAAAAAAAA");
  data.emplace_back("TTTTTTTTTTTTTTT");
  data.emplace_back("AAAAAAAAAAAAAAA");
  auto compressed = shared_tree::create_balanced(data);

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

  if (!errors) std::cout << "<Tree transposition> Finished without errors\n";
  return errors;
}

auto test_tree_iteration() -> int {
  auto path = "data/chmpxx";
  auto data = read_genome(path);
  auto compressed = shared_tree::create_balanced(data);
  auto errors = 0;

  auto i = 0;
  for (const auto c : compressed) {
    if (compressed[i] != c) {
      std::cerr << "<Tree iteration> Test failed: compressed[i] != compressed[i] for i = " << i << " out of " << data.size() - 1 << '\n'
        << "\tcompressed[i]\t= " << compressed[i] << '\n'
        << "\tcompressed[i]\t= " << c << '\n';
      ++errors;
    }
    ++i;
  }

  if (!errors) std::cout << "<Tree iteration> Finished without errors\n";
  return errors;
}

int main(int argc, char* argv[]) {
  auto errors = test_pointer() + test_chunks() + test_file_reader() + test_node() + test_tree_factory() + test_tree_similarity_transforms() + test_tree_iteration();
  if (errors) std::cerr << "Not all tests passed\n";
  return errors;
}