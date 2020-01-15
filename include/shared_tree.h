/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 */

#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <optional>
#include <stack>
#include <unordered_map>
#include <utility>
#include <vector>

#include "robin_hood.h"

#include "dna.h"
#include "utility.h"

/**
 *  Forward references.
 */
class node;
class pointer;

/******************************************************************************
 *  Annotated pointer to another node in the tree.
 *  Instead of references to leaf nodes, their values are stored directly.
 */
class pointer {
public:
  pointer(std::size_t index, bool mirror, bool transpose)
  : information{index}, mirror_bit{mirror}
  , transpose_bit{transpose}, leaf_bit{false} {}

  pointer(dna dna, bool mirror, bool transpose)
  : information{dna.to_ullong()}, mirror_bit{mirror}
  , transpose_bit{transpose}, leaf_bit{true} {}

  pointer(std::nullptr_t) : information{0}, mirror_bit{false}, transpose_bit{false}, leaf_bit{false} {}

  bool operator==(const pointer& other) const noexcept { return to_ullong() == other.to_ullong(); }
  bool operator!=(const pointer& other) const noexcept { return to_ullong() != other.to_ullong(); }
  bool operator<(const pointer& other) const noexcept { return information < other.information; }

  auto is_leaf() const noexcept -> bool { return leaf_bit; }
  auto is_mirrored() const noexcept -> bool { return mirror_bit; }
  auto is_transposed() const noexcept -> bool { return transpose_bit; }

  auto empty() const noexcept -> bool { return *this == nullptr; }
  auto data() const -> std::uint64_t { return information; }
  auto to_ullong() const noexcept -> unsigned long long;
  auto leaf() const -> dna;
  auto index() const -> std::uint64_t;

  auto mirrored() const noexcept { auto result = *this; result.mirror_bit ^= true; return result; }
  auto transposed() const noexcept { auto result = *this; result.transpose_bit ^= true; return result; }

private:
  std::uint64_t information : 60;
  bool mirror_bit : 1;
  bool transpose_bit : 1;
  bool leaf_bit : 1;
};

static inline auto operator<<(std::ostream& os, const pointer& p) -> std::ostream& {
  if (p.empty()) return os << "empty";
  if (p.is_leaf()) return os << "leaf (" << p.is_mirrored() << p.is_transposed() << "): " << p.leaf();
  else return os << "node (" << p.is_mirrored() << p.is_transposed() << "): " << p.index();
}


/******************************************************************************
 *  Node in a shared tree. Children are either other nodes or leaf nodes,
 *  containing the actual data stored.
 */
class node {
public:
  node(pointer left, pointer right) : children{left, right} {}
  node(pointer left) : children{left, nullptr} {}
  node(dna left, dna right) : children{left.canonical(), right.canonical()} {}
  node(dna left) : children{left.canonical(), nullptr} {}

  auto left() const -> const pointer& { return children[0]; }
  auto right() const -> const pointer& { return children[1]; }
  auto empty() const noexcept -> std::size_t { return children[0].empty() && children[1].empty(); }

  bool operator==(const node& other) const noexcept;
  bool operator!=(const node& other) const noexcept { return !(*this == other); }
  auto transformations(const node& other) const noexcept -> std::pair<bool, bool>;

  bool verbose_compare(const node& other) const noexcept;

private:
  std::array<pointer, 2> children;
};

static inline auto operator<<(std::ostream& os, const node& n) -> std::ostream& {
  return os << "node<" << n.left() << ", " << n.right() << ">";
}


/******************************************************************************
 *  Hash functions for the node class, so that it can be used in hash tables.
 */
namespace std {
  template<> struct hash<node> {
    auto operator()(const node& n) const noexcept -> std::size_t {
      const auto left = (n.left().to_ullong() & 0xfffffffffffffff) | static_cast<unsigned long long>(n.left().is_leaf()) << 62;
      const auto right = (n.right().to_ullong() & 0xfffffffffffffff) | static_cast<unsigned long long>(n.right().is_leaf()) << 62;
      const auto transposed = n.left().is_transposed() ^ n.right().is_transposed();
      const auto mirrored = n.left().is_mirrored() ^ n.left().is_mirrored();
      if (left < right) return detail::hash(1, transposed, mirrored, left, right);
      else return detail::hash(0, transposed, mirrored, right, left);
    }
  };
}

/******************************************************************************
 *  Genome, constructed as a balanced shared tree.
 */
class shared_tree {
public:
  shared_tree() : root{nullptr} {}
  shared_tree(const shared_tree&) = delete; // The self-referencing inside unordered set prevents copies
  shared_tree(shared_tree&& other) = default;

  template<typename Iterable>
  static auto create_balanced(Iterable&& data) -> shared_tree;

  auto node_count() const -> std::size_t { return nodes.size(); }
  auto width() const -> std::size_t { return children(root); }
  auto children(pointer parent) const -> std::size_t;
  void histogram(std::filesystem::path path) const;

  auto operator[](std::size_t index) const -> dna;
  auto access(pointer pointer) const -> const node&;

  struct iterator {
    iterator(const std::vector<node>& nodes, pointer root);

    auto operator*() noexcept -> dna;
    auto operator++() -> iterator&;
    auto operator!=(const iterator& other) const { return !stack.empty(); }

    auto access(pointer pointer) const -> const node&;
    void next_leaf();

    const std::vector<node>& nodes;
    std::vector<pointer> stack;
    std::vector<bool> transposed;
    bool mirrored = false;
  };

  using const_iterator = iterator;

  auto subtree(pointer parent) const {
    return iterator_pair(
      iterator{nodes, parent},
      iterator{nodes, nullptr}
    );
  }

  auto begin() { return iterator{nodes, root}; }
  auto end() { return iterator{nodes, nullptr}; }
  auto cbegin() { return const_iterator{nodes, root}; }
  auto cend() { return const_iterator{nodes, nullptr}; }

  template<typename... Args>
  void emplace_node(std::vector<pointer>& layer, Args&&... args);

private:
  template<typename Iterable>
  auto reduce_once(Iterable&& layer) -> std::vector<pointer>;

  template<typename Iterable>
  auto reduce(Iterable&& layer) -> pointer;

  pointer root;
  robin_hood::unordered_flat_map<node, std::size_t> indices;
  std::vector<node> nodes;
};

/**
 *  Constructs and emplaces a node inside the tree during its construction.
 *  In addition, the node is stored inside an auxiliary vector that helps
 *  store the structure of the layer being constructed, so that upper layers
 *  can also be constructed in full.
 */
template<typename... Args>
void shared_tree::emplace_node(std::vector<pointer>& layer, Args&&... args) {
  auto created_node = node{args...};
  auto insertion = indices.emplace(created_node, nodes.size()+1);
  auto canonical_node = (*insertion.first).first;
  auto index = (*insertion.first).second;

  auto mirrored = false;
  auto transposed = false;

  if (insertion.second) nodes.emplace_back(created_node);
  else {
    auto transform = created_node.transformations(canonical_node);
    mirrored = transform.first;
    transposed = transform.second;
    // if (mirrored || transposed) std::cout << created_node << ' ' << canonical_node << '\n';
  }

  layer.emplace_back(index, mirrored, transposed);
}

/**
 *  Reduces a layer, constructing the next by repeatedly creating pointers
 *  to nodes consisting of each adjacent pair of lower-level pointers.
 *  Nodes created in the process are stored inside the shared tree.
 */
template<typename Iterable>
auto shared_tree::reduce_once(Iterable&& layer) -> std::vector<pointer> {
  std::vector<pointer> next_layer;
  next_layer.reserve(layer.size()/2 + layer.size()%2);

  // Iterates over each pair, applying emplace_node.
  // If a singular element remains, applies emplace_node to that, too.
  foreach_pair(layer,
    [&](const auto& left, const auto& right) { emplace_node(next_layer, left, right); },
    [&](const auto& last) { emplace_node(next_layer, last); }
  );
  return next_layer;
}

/**
 *  Fully reduces a segment to its tree representation by repeated reduction of
 *  layers until one pointer remains. This pointer points to the root of the
 *  the resulting subtree, and is returned to the caller.
 */
template<typename Iterable>
auto shared_tree::reduce(Iterable&& segment) -> pointer {
  auto layer = reduce_once(segment);
  while (layer.size() > 1)
    layer = reduce_once(layer);
  return layer.front();
}

/**
 *  Constructs a shared binary tree from a range of data, using spatial
 *  subdivision for common subtree merging.
 *  Reduces the data into nodes segment-by-segment, reducing each segment in
 *  full before reducing the next. This reduces the memory requirements to
 *  those required for a single segment, at a small speed cost.
 */
template<typename Iterable>
auto shared_tree::create_balanced(Iterable&& data) -> shared_tree {
  constexpr auto segment_size = (1u<<29);
  auto result = shared_tree{};
  auto subroots = std::vector<pointer>{}; // Roots of all reduced segments

  for (auto segment : chunks(data, segment_size)) {
    auto subroot = result.reduce(segment);
    subroots.emplace_back(subroot);
  }

  result.root = result.reduce(subroots);
  return result;
}