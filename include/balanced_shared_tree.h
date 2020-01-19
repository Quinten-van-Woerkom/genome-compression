/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 *  Specialised for balanced trees, to improve compression speed and ratio for
 *  those trees. This does mean that unbalanced trees are not supported.
 */

#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <filesystem>
#include <utility>
#include <vector>

#include "robin_hood.h"

#include "dna.h"
#include "fasta_reader.h"
#include "utility.h"

namespace detail {
/****************************************************************************
 * class pointer:
 *  Pointer type used inside the balanced shared tree. Stores not just an
 *  index, but also two transformation bits to allow similar nodes to be
 *  stored as a single canonical node referred to with annotated pointers
 *  indicating the transformation applied to obtain the correct node.
 */
class pointer {
public:
  static constexpr auto address_bits = std::array{4, 12, 28, 60};

  pointer(std::nullptr_t = nullptr);
  pointer(const pointer& other, bool mirror = false, bool transpose = false);
  pointer(std::size_t index, bool mirror = false, bool transpose = false);
  pointer(dna leaf, bool mirror = false, bool transpose = false);

  bool operator==(const pointer& other) const noexcept { return to_ullong() == other.to_ullong(); }
  bool operator!=(const pointer& other) const noexcept { return to_ullong() != other.to_ullong(); }
  operator bool() const noexcept { return *this != nullptr; }

  bool empty() const noexcept { return *this == nullptr; }
  auto canonical() const noexcept { return data | ((std::uint64_t)segment << 62); }
  auto index() const noexcept -> std::size_t;
  auto leaf() const noexcept -> dna;

  auto bytes() const noexcept -> std::size_t { return (4 + address_bits[segment])/8; }
  void serialize(std::ostream& os) const;
  static auto deserialize(std::istream& is) -> pointer;

  bool is_mirrored() const noexcept { return mirror; }
  bool is_transposed() const noexcept { return transpose; }
  bool is_inverted() const noexcept { return mirror && transpose; }

  auto mirrored() const noexcept { return pointer{*this, true, false}; }
  auto transposed() const noexcept { return pointer{*this, false, true}; }
  auto inverted() const noexcept { return pointer{*this, true, true}; }

  /**
   * Returns the unsigned integer equivalent of the data stored in this
   * pointer.
   */
  auto to_ullong() const noexcept -> unsigned long long {
    return data | ((std::uint64_t)mirror << 60) | ((std::uint64_t)transpose << 61) | ((std::uint64_t)segment << 62);
  }

private:
  std::uint64_t data : 60;
  bool mirror : 1;
  bool transpose : 1;
  std::size_t segment : 2;
};

/****************************************************************************
 * class node:
 *  Node type representing inner nodes in the tree.
 *  Consists of two pointers to nodes or leaves one level down in the tree.
 */
class node {
public:
  node(pointer left, pointer right = nullptr) : children{left, right} {};

  auto left() const noexcept { return children[0]; }
  auto right() const noexcept { return children[1]; }

  auto mirrored() const noexcept { return node{children[1].mirrored(), children[0].mirrored()}; }
  auto transposed() const noexcept { return node{children[0].transposed(), children[1].transposed()}; }
  auto inverted() const noexcept { return mirrored().transposed(); }

  bool operator==(const node& other) const noexcept;
  bool operator!=(const node& other) const noexcept { return !(*this == other); };
  auto transformations(const node& other) const noexcept -> std::pair<bool, bool>;
  
  auto bytes() const noexcept { return left().bytes() + right().bytes(); }
  void serialize(std::ostream& os) const;
  static auto deserialize(std::istream& is) -> node;

private:
  std::array<pointer, 2> children;
};
}

/******************************************************************************
 * ostream operators
 *  Ostream output formatting for pointers and nodes.
 */
inline auto& operator<<(std::ostream& os, const detail::pointer& pointer) {
  if (pointer.empty()) return os << "empty";
  else return os << "index " << pointer.index();
}

inline auto& operator<<(std::ostream& os, const detail::node& node) {
  return os << "node<" << node.left() << ", " << node.right() << ">";
}


/******************************************************************************
 * std::hash<> specialisation
 *  Required for usage of nodes in a hashtable.
 */
namespace std {
  template<> struct hash<detail::node> {
    auto operator()(const detail::node& n) const noexcept -> std::size_t {
      const auto left = n.left().canonical();
      const auto right = n.right().canonical();
      const auto transposed = n.left().is_transposed() ^ n.right().is_transposed();
      const auto mirrored = n.left().is_mirrored() ^ n.right().is_mirrored();
      if (left < right) return detail::hash(1, transposed, mirrored, left, right);
      else return detail::hash(0, transposed, mirrored, right, left);
    }
  };
}

/******************************************************************************
 * class balanced_shared_tree:
 *  Shared binary tree class that exploits structural properties of balanced
 *  trees to store its directed acyclic graph representation more effectively.
 */
class balanced_shared_tree {
public:
  using pointer = detail::pointer;
  using node = detail::node;

  balanced_shared_tree() = default;

  balanced_shared_tree(std::filesystem::path path)
  : balanced_shared_tree{fasta_reader{path}} {};

  balanced_shared_tree(fasta_reader file);
  balanced_shared_tree(std::vector<dna>& data);

  auto depth() const { return nodes.size(); }
  auto width() const { assert(nodes.back().size() == 1); return children(nodes.size()-1, root); }

  auto children(std::size_t layer, pointer pointer) const -> std::size_t;
  auto node_count() const -> std::size_t;
  auto node_count(std::size_t layer) const { return nodes[layer].size(); }

  auto access_leaf(std::size_t index, pointer pointer) const -> dna;
  auto access_node(std::size_t layer, pointer pointer) const -> node;
  auto operator[](std::uint64_t index) const -> dna;

  void add_layer() { nodes.emplace_back(); }
  void emplace_node(std::size_t layer, node node);
  void emplace_leaf(dna leaf);

  void histogram(std::filesystem::path) const;
  void print_unique(std::ostream& os) const;

  auto bytes() const noexcept -> std::size_t;
  void serialize(std::ostream& os) const;
  static auto deserialize(std::istream& is) -> balanced_shared_tree;

  struct iterator {
    using pointer = balanced_shared_tree::pointer;
    using node = balanced_shared_tree::node;

    struct status {
      status(std::size_t layer, pointer current)
      : layer{layer}, current{current} {}

      std::size_t layer;  // Leaf node is denoted by max std::size_t value
      iterator::pointer current;
    };

    iterator(balanced_shared_tree& nodes, std::size_t layer, pointer root);

    auto operator*() const noexcept -> dna;
    auto operator++() -> iterator&;
    auto operator!=(const iterator& other) const { return !stack.empty(); }
    void next_leaf();

    balanced_shared_tree& parent;
    std::vector<status> stack;
  };

  using const_iterator = iterator;

  auto begin() { return iterator{*this, nodes.size()-1, root}; }
  auto end() { return iterator{*this, 0, nullptr}; }

private:
  std::vector<std::vector<node>> nodes;
  pointer root;
};


/******************************************************************************
 * class tree_constructor:
 *  Helper class in construction of a balanced shared tree.
 *  Contains the maps used to link nodes to pointers or leaves.
 */
class tree_constructor {
public:
  using pointer = balanced_shared_tree::pointer;
  using node = balanced_shared_tree::node;
  using hash_map = robin_hood::unordered_flat_map<node, std::uint64_t>;

  tree_constructor(balanced_shared_tree& parent);

  void emplace(std::vector<pointer>& layer, std::size_t layer_index, pointer left, pointer right = nullptr);

  template<typename Iterable>
  auto reduce_layer(Iterable&& layer, std::size_t index) -> std::vector<pointer>;

  template<typename Iterable>
  auto reduce_segment(Iterable&& segment, std::size_t index = 0) -> pointer;

  template<typename Iterable>
  auto reduce(Iterable&& data) -> pointer;

private:
  balanced_shared_tree& parent;
  std::vector<hash_map> nodes;
  std::vector<pointer> roots;
};

/**
 * Reduces the current layer by constructing nodes out of the pointers it
 * contains, and emplacing those nodes in the correct layers and maps, if
 * necessary.
 */
template<typename Iterable>
auto tree_constructor::reduce_layer(Iterable&& iterable, std::size_t index) -> std::vector<pointer> {
  auto layer = std::vector<pointer>{};
  layer.reserve(iterable.size()/2 + iterable.size()%2);
  if (parent.depth() <= index) {
    parent.add_layer();
    nodes.emplace_back();
  }

  foreach_pair(iterable,
    [&](auto left, auto right) { emplace(layer, index, left, right); },
    [&](auto last) { emplace(layer, index, last); }
  );

  return layer;
}

/**
 * Reduces the layer, saving any encountered nodes and leaves in the shared
 * tree according to canonical representation.
 */
template<typename Iterable>
auto tree_constructor::reduce_segment(Iterable&& segment, std::size_t index) -> pointer {
  auto layer = reduce_layer(segment, index++);
  for (; layer.size() > 1; ++index)
    layer = reduce_layer(layer, index);
  return layer.front();
}

/**
 * Reduces the data by dividing it into segments, each of which is fully
 * reduced. The resulting tree roots are then also reduced to obtain the
 * final tree representation.
 */
template<typename Iterable>
auto tree_constructor::reduce(Iterable&& data) -> pointer {
  constexpr auto segment_size = (1u<<10);
  auto segment_roots = std::vector<pointer>{};

  for (auto segment : chunks(data, segment_size)) {
    auto segment_root = reduce_segment(segment);
    segment_roots.emplace_back(segment_root);
  }
  return reduce_segment(segment_roots, nodes.size());
}