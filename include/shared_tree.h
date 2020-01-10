/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 */

#pragma once

#include <cassert>
#include <memory>
#include <memory_resource>
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

  auto is_leaf() const noexcept -> bool { return leaf == true; }
  auto empty() const noexcept -> bool { return !is_leaf() && subnode == nullptr; }
  auto size() const noexcept -> std::size_t;
  
  auto get_leaf() const -> const dna&;
  auto get_node() const -> const node&;

private:
  bool leaf : 1;
  
  union {
    const void* void_pointer;
    const node* subnode;
    dna data;
  };
};

static auto operator<<(std::ostream& os, const pointer& p) -> std::ostream& {
    if (p.empty()) return os << "empty";
    else if (p.is_leaf()) return os << "leaf: " << p.get_leaf();
    else return os << "node: " << &p.get_node();
  }

/******************************************************************************
 *  Node in a shared tree. Children are either other nodes or leaf nodes,
 *  containing the actual data stored.
 */
class node {
public:
  using pointer = pointer;

  node(pointer left, pointer right) : left_{left}, right_{right} {}
  node(pointer left) : left_{left}, right_{nullptr} {}
  
  auto left() const -> const pointer& { return left_; }
  auto right() const -> const pointer& { return right_; }
  
  auto hash() const noexcept -> std::size_t { return detail::hash(left_, right_); }
  auto size() const noexcept -> std::size_t { return left_.size() + right_.size(); }
  auto empty() const noexcept -> bool { return left_.empty() && right_.empty(); }

  bool operator==(const node& other) const noexcept {
    return left_ == other.left_ && right_ == other.right_;
  }

private:
  pointer left_, right_;
};

static auto operator<<(std::ostream& os, const node& n) -> std::ostream& {
  return os << "node<" << n.left() << ", " << n.right() << ">";
}

namespace std {
  template<> struct hash<node> {
    auto operator()(const node& n) const noexcept -> std::size_t {
      return n.hash();
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

  // Number of nodes stored
  auto node_count() const -> std::size_t { return nodes.size(); }
  auto width() const -> std::size_t { return root.size(); }

  auto operator[](std::size_t index) const -> const dna&;
  
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
 *  Constructs a shared binary tree from a range of data, using spatial
 *  subdivision for common subtree merging.
 *  Reduces the data into nodes layer-by-layer, reducing each layer in segments
 *  that fit in memory.
 */
template<typename Iterable>
auto shared_tree::create_balanced(Iterable&& data) -> shared_tree {
  constexpr auto segment_size = (1u<<29);
  auto result = shared_tree{};
  auto segments = std::vector<pointer>{};

  auto reduce_once = [&](auto&& layer) {
    std::vector<pointer> next_layer;
    next_layer.reserve(layer.size()/2 + layer.size()%2);

    // Emplaces a node in the next layer, using arguments passed
    auto emplace_node = [&](const auto&... args) {
      auto created_node = node{args...};
      auto insertion = result.nodes.emplace(created_node);
      auto& canonical_node = *(insertion.first);
      next_layer.emplace_back(canonical_node);
    };

    // Iterates over each pair, applying emplace_node.
    // If a singular element remains, applies emplace_node to that, too.
    foreach_pair(layer, emplace_node, emplace_node);
    return next_layer;
  };

  // Reduces a layer repeatedly, returning the root node of the resulting tree
  auto reduce = [&](auto&& layer) {
    while (layer.size() > 1)
      layer = reduce_once(layer);
    return layer.front();
  };

  for (auto segment : chunks(data, segment_size)) {
    auto layer = reduce_once(segment);
    auto segment_root = reduce(layer);
    segments.emplace_back(segment_root);
  }

  result.root = reduce(segments);
  return result;
}