/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 */

#pragma once

#include <cassert>
#include <optional>
#include <stack>
#include <unordered_set>
#include <utility>
#include <vector>

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
  pointer(std::nullptr_t = nullptr) : leaf{false}, void_pointer{nullptr} {}  // Special construct to denote the absence of data
  pointer(const node& subnode) : leaf{false}, subnode{&subnode} {}
  pointer(dna data) : leaf{true}, data{data} {}

  operator std::size_t() const noexcept { return detail::hash(reinterpret_cast<std::size_t>(void_pointer), leaf); }

  bool operator==(const pointer& other) const noexcept {
    return leaf == other.leaf && void_pointer == other.void_pointer;
  }

  // Accessor
  auto operator[](std::size_t index) const -> const dna&;

  auto is_leaf() const noexcept -> bool { return leaf == true; }
  auto empty() const noexcept -> bool { return !is_leaf() && subnode == nullptr; }
  auto size() const noexcept -> std::size_t;
  
  auto get_leaf() const -> const dna&;
  auto get_node() const -> const node&;

  friend auto operator<<(std::ostream& os, const pointer& p) -> std::ostream& {
    if (p.empty()) return os << "empty";
    else if (p.is_leaf()) return os << "leaf: " << p.data;
    else return os << "node: " << p.subnode;
  }

private:
  bool leaf : 1;
  
  union {
    const void* void_pointer;
    const node* subnode;
    dna data;
  };
};

/******************************************************************************
 *  Node in a non-owning tree. Children are either other nodes or leaf nodes,
 *  containing the actual data stored.
 */
class node {
public:
  using pointer = pointer;

  node(pointer left, pointer right) : left_{left}, right_{right} {}
  node(pointer left) : left_{left}, right_{nullptr} {}
  
  auto left() const -> pointer { return left_; }
  auto right() const -> pointer { return right_; }
  
  auto hash() const noexcept -> std::size_t { return detail::hash(left_, right_); }
  auto size() const noexcept -> std::size_t { return left_.size() + right_.size(); }
  auto full() const noexcept -> bool { return !left_.empty() && !right_.empty(); }
  auto empty() const noexcept -> bool { return left_.empty() && right_.empty(); }

  // Accessor
  auto operator[](std::size_t index) const -> const dna& {
    if (index < left_.size()) return left_[index];
    else return right_[index - left_.size()];
  }

  bool operator==(const node& other) const noexcept {
    return left_ == other.left_ && right_ == other.right_;
  }

  friend auto operator<<(std::ostream& os, const node& n) -> std::ostream& {
    return os << "node<" << n.left_ << ", " << n.right_ << ">";
  }

private:
  pointer left_, right_;
};

namespace std {
  template<> struct hash<node> {
    auto operator()(const node& n) const noexcept -> std::size_t {
      return n.hash();
    }
  };
}

/**
 *  Genome, constructed as a balanced shared tree.
 *  Note that the tree itself is non-owning; rather, the nodes must be stored
 *  elsewhere to allow for canonicalisation of nodes.
 *  Precondition: each node is assumed to exist only once.
 */
class shared_tree {
public:
  shared_tree() : root{nullptr} {}
  shared_tree(const shared_tree&) = delete; // Self-referencing inside unordered set prevents copies
  shared_tree(shared_tree&&) = default;

  template<typename Iterable>
  static auto create_balanced(Iterable&& data) -> shared_tree;

  // Number of nodes stored
  auto size() const -> std::size_t { return nodes.size(); }

  // Length of the data sequence
  auto width() const -> std::size_t { return root.size(); }

  auto operator[](std::size_t index) const -> const dna& {
    return root[index];
  }
  
  void print_unique() const {
    for (const auto& node : nodes)
        std::cout << node << '\n';
  }

  class iterator {
  public:
    iterator(shared_tree& parent, std::size_t index = 0) : parent{parent}, index{index} {}
    auto& operator*() { return parent[index]; }
    auto& operator++() { ++index; return *this; }
    auto operator++(int) { auto temp = *this; ++index; return temp; }
    auto& operator--() { --index; return *this; }
    auto operator--(int) { auto temp = *this; --index; return temp; }
    auto operator!=(const iterator& other) const { return &parent != &other.parent || index != other.index; }

  private:
    shared_tree& parent;
    std::size_t index;
  };

  using const_iterator = iterator;

  auto begin() { return iterator{*this, 0}; }
  auto end() { return iterator{*this, root.size()}; }
  auto cbegin() { return const_iterator{*this, 0}; }
  auto cend() { return const_iterator{*this, root.size()}; }

private:
  pointer root;
  std::unordered_set<node> nodes;
};

/**
 *  Constructs a shared binary tree from a vector of data, using spatial
 *  subdivision for common subtree merging.
 *  Reduces the data into nodes layer-by-layer, reducing each layer in segments
 *  that fit in memory.
 */
template<typename Iterable>
auto shared_tree::create_balanced(Iterable&& data) -> shared_tree {
  constexpr auto segment_size = (1u<<26);
  auto result = shared_tree{};
  auto segments = std::vector<pointer>{};

  auto reduce_layer = [&](const auto& layer) {
    std::vector<pointer> next_layer;
    next_layer.reserve(layer.size()/2 + layer.size()%2);
    pairwise_apply(
      layer,
      [&](const auto& left, const auto& right) {
        auto created_node = node{left, right};
        auto insertion = result.nodes.emplace(created_node);
        auto& canonical_node = *(insertion.first);
        next_layer.emplace_back(canonical_node);
      },
      [&](const auto& left) {
        auto created_node = node{left};
        auto insertion = result.nodes.emplace(created_node);
        auto& canonical_node = *(insertion.first);
        next_layer.emplace_back(canonical_node);
      }
    );
    return next_layer;
  };

  for (const auto&& segment : chunks(data, segment_size)) {
    auto layer = reduce_layer(segment);

    while (layer.size() > 1) {
      auto next_layer = reduce_layer(layer);
      std::swap(layer, next_layer);
    }
    segments.emplace_back(layer.front());
  }

  while (segments.size() > 1) {
    auto next_layer = reduce_layer(segments);
    std::swap(next_layer, segments);
  }

  result.root = pointer{segments.front()};
  return result;
}