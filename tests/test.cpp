/**
 *  Unit tests for the implementation.
 */

#include <array>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "balanced_shared_tree.h"
#include "dna.h"
#include "fasta_reader.h"
#include "shared_tree.h"
#include "utility.h"

auto test_pointer() -> int {
  auto errors = 0;
  auto leaf = dna{"TTGAACGAGAAGCCG"};
  auto basis = detail::pointer{leaf};
  auto transposed = basis.transposed();
  auto mirrored = basis.mirrored();
  auto inner_pointer = detail::pointer{3280, true, true};

  if (basis.leaf() != leaf) {
    std::cerr << "<Tree pointer> DNA conversion does not match input\n";
    ++errors;
  }

  if (basis == transposed) {
    std::cerr << "<Tree pointer> Non-null pointers can never match when transposed\n";
    ++errors;
  }

  if (basis == mirrored) {
    std::cerr << "<Tree pointer> Pointers that are not invariant under mirroring should not match when mirrored\n";
    ++errors;
  }

  if (inner_pointer.index() != 3280) {
    std::cerr << "<Tree pointer> Index conversion error: "
      << inner_pointer.index() << " != " << 3280 << '\n';
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
  // auto path = "../GRCm38.fna";
  auto data = read_genome(path);
  // auto compressed = shared_tree::create_balanced(data);
  auto compressed = balanced_shared_tree{path};
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
  auto data = dna{"CCTCTGCCTCTGCCT"};
  auto other = dna{"CTGCCTCTGCCTCTG"};
  auto basis = data;
  auto transposed = data.transposed();
  auto mirrored = data.mirrored();
  auto both = data.transposed().mirrored();
  auto con_layer = std::vector<dna>{};
  auto layer = std::vector<detail::pointer>{};

  con_layer.emplace_back(basis);
  con_layer.emplace_back(transposed);
  con_layer.emplace_back(both);
  con_layer.emplace_back(mirrored);
  auto tree = balanced_shared_tree{con_layer};

  if (node{basis, transposed} != node{both, mirrored}
  || std::hash<node>()(node{basis, transposed}) != std::hash<node>()(node{both, mirrored})) {
    std::cerr << "<Tree similarity> Mirrored nodes should compare to be equal\n";
  }

  auto basis_node = node{basis, basis};
  auto mirrored_node = node{mirrored, mirrored};
  auto transposed_node = node{transposed, transposed};
  auto mix_node = node{basis, transposed};
  auto mirrored_mix_node = node{both, mirrored};

  if (basis_node != mirrored_node
    || std::hash<node>()(basis_node) != std::hash<node>()(mirrored_node)) {
    std::cerr << "<Tree similarity> Mirrored nodes should compare to be equal\n";
    ++errors;
  }

  if (basis_node != transposed_node
    || std::hash<node>()(basis_node) != std::hash<node>()(transposed_node)) {
    std::cerr << "<Tree similarity> Transposed nodes should compare to be equal\n";
    ++errors;
  }

  if (mix_node != mirrored_mix_node
    || std::hash<node>()(mix_node) != std::hash<node>()(mirrored_mix_node)) {
    std::cerr << "<Tree similarity> Transposed and mirrored nodes should compare to be equal\n";
    ++errors;
  }

  auto left = node{basis, other};
  auto right = node{other, basis};
  if (left == right || std::hash<node>()(left) == std::hash<node>()(right)) {
    std::cerr << "<Tree similarity> Nodes with reversed children that are not mirrored should not match if they are not invariant under transformation\n";
    ++errors;
  }

  if (tree.node_count() != 3) {
    std::cerr << "<Tree similarity> Similar nodes do not merge: node count is "
      << tree.node_count() << ", not matching the expected 3:\n";
    tree.print_unique(std::cerr);
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
  auto compressed = balanced_shared_tree(data);

  if (data.size() != compressed.width()) {
    std::cerr << "<Tree transposition> Test failed: raw data size (" << data.size()
      << ") does not match compressed data size (" << compressed.width() << ", " << compressed.node_count() << " nodes)\n";
    ++errors;
  }

  auto i = 0;
  for (const auto c : compressed) {
    if (data[i] != c) {
      std::cerr << "<Tree transposition> Test failed: data[i] != compressed[i] for i = " << i << " out of " << data.size() - 1 << '\n'
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
  // auto path = "../GRCm38.fna";
  auto data = read_genome(path);
  auto compressed = balanced_shared_tree{path};
  auto errors = 0;

  auto i = 0;
  for (const auto c : compressed) {
    if (data[i] != c) {
      std::cerr << "<Tree iteration> Test failed: data[i] != compressed[i] for i = " << i << " out of " << data.size() - 1 << '\n'
        << "\tdata[i]\t\t= " << data[i] << '\n'
        << "\tcompressed[i]\t= " << c << '\n';
      ++errors;
    }
    ++i;
  }

  i = 0;
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
  auto errors = test_pointer() + test_chunks() + test_file_reader() + test_node() + test_tree_factory() + test_tree_similarity_transforms() + test_tree_transposition() + test_tree_iteration();
  if (errors) std::cerr << "Not all tests passed\n";
  return errors;
}