/**
 *  Unit tests for the implementation.
 */

#include <array>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include "balanced_shared_tree.h"
#include "dna.h"
#include "fasta_reader.h"
#include "utility.h"

using detail::node;
using detail::pointer;

#define TEST_START(name) \
    auto errors = 0;\
    auto expects = [&](bool test, auto... message) {\
        if (!test) {\
        ++errors;\
        std::cout << "<" << name << "> ";\
        ((std::cout << message), ...);\
        std::cout << '\n';\
        }\
    };

#define TEST_END(name) \
    if (!errors) std::cout << "<" << name << "> Finished without errors\n"; \
    return errors;

auto test_pointer() -> int {
  TEST_START("Tree pointer");

  auto leaf = dna::random(3);
  auto basis = detail::pointer{leaf};
  auto transposed = basis.transposed();
  auto mirrored = basis.mirrored();
  auto inner_pointer = detail::pointer{3280, true, true};

  expects(basis.leaf() == leaf, "DNA conversion does not match input");
  expects(basis != transposed, "Non-null pointers can never match when transposed");
  expects(basis != mirrored, "Pointers that are not invariant under mirroring should not match when mirrored");
  expects(inner_pointer.index() == 3280, "Index conversion error: ", inner_pointer.index(), " != ", 3280);

  TEST_END("Tree pointer");
}

auto test_chunks() -> int {
  TEST_START("Chunks");

  auto a = std::array{1, 2, 3, 4, 5, 6, 7, 8};
  auto index = 1;

  for (auto chunk : chunks(a, 2)) {
    for (auto v : chunk) {
      expects(v == index, "Value mismatch: ", v, " != ", index);
      ++index;
    }
  }
  
  TEST_END("Chunks");
}

auto test_file_reader() -> int {
  TEST_START("File reader");

  auto path = "data/chmpxx";
  auto size = std::filesystem::file_size(path);
  auto buffered = fasta_reader{path};
  auto direct = std::vector<std::array<char, dna::size()>>(size/dna::size());
  auto file = std::ifstream{path, std::ios::binary};

  file.read(reinterpret_cast<char*>(direct.data()), size);

  expects(
    direct.size() == size/dna::size(),
    "File path size does not match read size within accuracy bounds\n",
    "File path: ", size, '\n',
    "Direct read: ", direct.size()*dna::size()
  );

  auto i = 0u;
  for (auto b : buffered) {
    auto d = std::string_view{direct[i].data(), direct[i].size()};
    expects(
      b == d,
      "buffered[i] != direct[i] for i = ", i, " out of ", direct.size() - 1, '\n',
      "\tbuffered[i]\t= ", b, '\n',
      "\tdirect[i]\t= ", d
    );
    ++i;
  }

  expects(
    i == size/dna::size(),
    "File path size does not match buffered read size within accuracy bounds\n",
    "File path: ", size, '\n',
    "Buffered read: ", i*dna::size());

  TEST_END("File reader");
}

auto test_tree_factory() -> int {
  TEST_START("Tree factory");

  auto path = "data/chmpxx";
  // auto path = "../GRCm38.fna";
  auto data = read_genome(path);
  auto compressed = balanced_shared_tree{path};

  expects(
    data.size() == compressed.width(),
    "Raw data size (", data.size(), ") does not match compressed data size (",
    compressed.width(), ", ", compressed.node_count(), "nodes)"
  );

  auto i = 0;
  for (const auto c : compressed) {
    expects(
      data[i] == c,
      "data[i] != compressed[i] for i = ", i, " out of ", data.size() - 1, '\n',
      "\tdata[i]\t\t= ", data[i], '\n',
      "\tcompressed[i]\t= ", c
    );
    ++i;
  }

  TEST_END("Tree factory");
}

auto test_similarity_transforms() -> int {
  TEST_START("Similarity transforms");

  {
    auto basis = dna::random(1);
    auto transposed = basis.transposed();
    auto mirrored = basis.mirrored();
    auto inverted = basis.inverted();

    auto basis_node = node{basis, basis};
    auto mirrored_node = node{mirrored, mirrored};
    auto transposed_node = node{transposed, transposed};
    auto mix_node = node{basis, transposed};
    auto mirrored_mix_node = node{inverted, mirrored};

    expects(basis_node == mirrored_node, "Mirrored nodes should compare to be equal");
    expects(basis_node == transposed_node, "Transposed nodes should compare to be equal");
    expects(mix_node == mirrored_mix_node, "Inverted nodes should compare to be equal");
    expects(std::hash<node>()(basis_node) == std::hash<node>()(mirrored_node), "Mirrored nodes should have the same hash");
    expects(std::hash<node>()(basis_node) == std::hash<node>()(transposed_node), "Transposed nodes should have the same hash");
    expects(std::hash<node>()(mix_node) == std::hash<node>()(mirrored_mix_node), "Inverted nodes should have the same hash");
  }

  {
    auto basis = dna::random(1);
    auto other = dna::random(42);
    auto left = node{basis, other};
    auto right = node{other, basis};

    expects(left != right, "Nodes with reversed children that are not mirrored should not match if they are not invariant under transformation");
  }

  {
    auto basis = dna::random(1);
    auto transposed = basis.transposed();
    auto mirrored = basis.mirrored();
    auto inverted = basis.inverted();

    auto con_layer = std::vector{basis, transposed, inverted, mirrored};
    auto tree = balanced_shared_tree{con_layer};

    expects(
        tree.node_count() == 3,
        "Similar modes should merge: node count is ", tree.node_count(),
        ", not matching the expected 3"
    );
  }

  TEST_END("Similarity transforms");
}

auto test_tree_transposition() -> int {
  TEST_START("Tree transposition");

  auto a = dna::random(0);
  auto t = a.transposed();
  auto data = std::vector<dna>{a, a, t, a};
  auto compressed = balanced_shared_tree(data);

  expects(
    data.size() == compressed.width(),
    "Raw data size (", data.size(), ") does not match compressed data size (",
    compressed.width(), ", ", compressed.node_count(), " nodes)"
  );

  auto i = 0;
  for (const auto c : compressed) {
    expects(
      data[i] == c,
      "data[i] != compressed[i] for i = ", i, " out of ", data.size() - 1, '\n',
      "\tdata[i]\t\t= ", data[i], '\n',
      "\tcompressed[i]\t= ", c
    );
    ++i;
  }

  TEST_END("Tree transposition");
}

auto test_serialization() -> int {
  TEST_START("Serialization");

  {
    auto a = dna::random(0);
    auto b = dna::random(1);
    auto c = dna::random(2);

    auto basis = pointer{a};
    auto other = pointer{b};
    auto left = node{basis, other};
    auto right = node{other, basis};
    auto stream = std::stringstream{};

    basis.serialize(stream);
    other.serialize(stream);
    left.serialize(stream);
    right.serialize(stream);

    auto sbasis = pointer::deserialize(stream);
    auto sother = pointer::deserialize(stream);
    auto sleft = node::deserialize(stream);
    auto sright = node::deserialize(stream);

    expects(basis == sbasis, "Serialization and deserialization should result in identical pointers: ", basis, " != ", sbasis);
    expects(other == sother, "Serialization and deserialization should result in identical pointers: ", other, " != ", sother);
    expects(left == sleft, "Serialization and deserialization should result in identical nodes: ", left, " != ", sleft);
    expects(right == sright, "Serialization and deserialization should result in identical nodes: ", right, " != ", sright);

    auto data = std::vector{a, b, b, a, c, a};
    auto tree = balanced_shared_tree{data};
    stream = std::stringstream{};
    tree.serialize(stream);
    auto load = balanced_shared_tree::deserialize(stream);

    expects(tree.width() == load.width(), "Serialization and deserialization should result in identical tree size");

    for (auto i = 0u; i < tree.width(); ++i) {
      expects(tree[i] == load[i], "Serialization and deserialization should result in identical tree: ", tree[i], " != ", load[i]);
    }
  }

  TEST_END("Serialization");
}

auto test_tree_iteration() -> int {
  TEST_START("Tree iteration");

  auto path = "data/chmpxx";
  auto data = read_genome(path);
  auto compressed = balanced_shared_tree{path};

  auto i = 0;
  for (const auto c : compressed) {
    expects(
      data[i] == c,
      "data[i] != compressed[i] for i = ", i, " out of ", data.size() - 1, '\n',
      "\tdata[i]\t\t= ", data[i], '\n',
      "\tcompressed[i]\t= ", c
    );
    ++i;
  }

  i = 0;
  for (const auto c : compressed) {
    expects(
      compressed[i] == c,
      "compressed[i] != compressed[i] for i = ", i, " out of ", data.size() - 1, '\n',
      "\tcompressed[i]\t\t= ", data[i], '\n',
      "\tcompressed[i]\t= ", c
    );
    ++i;
  }

  TEST_END("Tree iteration");
}

int main(int argc, char* argv[]) {
  auto errors = test_pointer() + test_chunks() + test_file_reader() + test_tree_factory() + test_similarity_transforms() + test_tree_transposition() + test_serialization() + test_tree_iteration();
  if (errors) std::cerr << "Not all tests passed\n";
  return errors;
}