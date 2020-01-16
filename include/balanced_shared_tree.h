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
#include <filesystem>
#include <vector>

#include "robin_hood.h"

#include "dna.h"
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
  pointer(std::uint64_t index, bool mirror, bool transpose)
  : data{index}, mirror{mirror}, transpose{transpose} {
    // Ensures that a transformed nullptr is also null
    const auto is_null = data == (1ull<<63)-1;
    mirror |= is_null;
    transpose |= is_null;
  }

  pointer(std::nullptr_t = nullptr) : pointer{(1ull<<63)-1, true, true} {}

  bool operator==(const pointer& other) const noexcept { return to_ullong() == other.to_ullong(); }
  bool operator!=(const pointer& other) const noexcept { return to_ullong() != other.to_ullong(); }
  operator bool() const noexcept { return *this != nullptr; }

  bool empty() const noexcept { return *this == nullptr; }
  auto index() const noexcept { assert(!empty()); return data; }

  bool is_mirrored() const noexcept { return mirror; }
  bool is_transposed() const noexcept { return transpose; }
  bool is_inverted() const noexcept { return mirror && transpose; }

  auto mirrored() const noexcept { return pointer{data, !mirror, transpose}; }
  auto transposed() const noexcept { return pointer{data, mirror, !transpose}; }
  auto inverted() const noexcept { return pointer{data, !mirror, !transpose}; }

  /**
   * Returns the unsigned integer equivalent of the data stored in this
   * pointer.
   */
  auto to_ullong() const noexcept -> unsigned long long {
    return data | ((std::uint64_t)mirror << 62) | ((std::uint64_t)transpose << 63);
  }

private:
  std::uint64_t data : 62;
  bool mirror : 1;
  bool transpose : 1;
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

  auto mirrored() const noexcept { return node{children[0].mirrored(), children[1].mirrored()}; }
  auto transposed() const noexcept { return node{children[0].transposed(), children[1].transposed()}; }
  auto inverted() const noexcept { return mirrored().transposed(); }

  bool operator==(const node& other) const noexcept;
  bool operator!=(const node& other) const noexcept { return !(*this == other); };
  auto transformations(const node& other) const noexcept -> std::pair<bool, bool>;

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
      const auto left = n.left().index();
      const auto right = n.left().index();
      const auto transposed = n.left().is_transposed() ^ n.right().is_transposed();
      const auto mirrored = n.left().is_mirrored() ^ n.left().is_mirrored();
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

  balanced_shared_tree(std::filesystem::path path);

  auto depth() const { return nodes.size(); }
  auto width() const { return children(nodes.size()-1, root); }

  auto children(std::size_t layer, pointer pointer) const -> std::size_t;
  auto node_count() const -> std::size_t;
  auto node_count(std::size_t layer) const { return nodes[layer].size(); }
  auto leaf_count() const { return leaves.size(); }

  auto access_leaf(pointer pointer) const -> dna;
  auto access_node(std::size_t layer, pointer pointer) const -> node;
  auto operator[](std::uint64_t index) const -> dna;

  void add_layer() { nodes.emplace_back(); }
  void emplace_node(std::size_t layer, node node);
  void emplace_leaf(dna leaf);

  void histogram(std::filesystem::path) {};

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
  std::vector<dna> leaves;
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
  
  // template<typename T>
  // using hash_map = robin_hood::unordered_flat_map<T, std::uint64_t>;
  template<typename T>
  using hash_map = robin_hood::unordered_node_map<T, std::uint64_t>;

  tree_constructor(balanced_shared_tree& parent);

  auto root() const -> pointer { return roots.front(); }

  void emplace_node(pointer left, pointer right = nullptr);
  void emplace_leaf(const dna& leaf);
  void reduce_once();

  template<typename Iterable>
  void reduce(Iterable&& segment);
  void reduce();

private:
  balanced_shared_tree& parent;
  std::vector<hash_map<node>> nodes;
  hash_map<dna> leaves;
  std::vector<pointer> current_layer;
  std::vector<pointer> next_layer;
  std::vector<pointer> roots;
  std::size_t layer = 0;
};

/**
 * Reduces the layer, saving any encountered nodes and leaves in the shared
 * tree according to canonical representation.
 */
template<typename Iterable>
void tree_constructor::reduce(Iterable&& segment) {
  layer = 0;
  for (auto leaf : segment) emplace_leaf(leaf);
  while (current_layer.size() > 1) reduce_once();
  roots.emplace_back(current_layer.front());
}