/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 */

#pragma once

#include <cassert>
#include <unordered_set>
#include <vector>

#include "utility.h"

/**
 *  Node in a non-owning tree. Children are either other nodes or leaf nodes,
 *  containing the actual data stored.
 */
class node {
public:
  class pointer;

  node(pointer left, pointer right) : left{left}, right{right} {}
  node(pointer left) : left{left}, right{nullptr} {}

  // TODO: Add checks for empty children
  auto left_leaf() const -> const dna& { assert(left.is_leaf()); return left.get_leaf(); }
  auto right_leaf() const -> const dna& { assert(right.is_leaf()); return right.get_leaf(); }
  auto left_node() const -> const node& { assert(!left.is_leaf()); return left.get_node(); }
  auto right_node() const -> const node& { assert(!right.is_leaf()); return right.get_node(); }
  
  // Hash function
  auto hash() const noexcept -> std::size_t { return detail::hash(left, right); }

  auto size() const noexcept -> std::size_t { return left.size() + right.size(); }

  // For now, only mirror symmetry exists.
  // In the future, more similarities might be added.
  auto similarities(const node& other) const {
    assert(*this == other);
    if (left == other.right && right == other.left) return true;
    else return false;
  }

  // Accessor
  auto operator[](std::size_t index) const -> const dna& {
    if (index < left.size()) return left[index];
    else return right[index - left.size()];
  }

  bool operator==(const node& other) const noexcept {
    return left == other.left && right == other.right;
  }

  friend auto operator<<(std::ostream& os, const node& n) -> std::ostream& {
    return os << "node<" << n.left << ", " << n.right << ">";
  }

  /**
   *  Annotated pointer to another node in the tree.
   *  This node can be either a leaf or an internal node.
   *  Data is assumed to be approximately the same size as a pointer, so that
   *  it makes sense to store the data instead of a pointer.
   */
  class pointer {
  public:
    pointer(std::nullptr_t = nullptr) : size_{0}, void_pointer{nullptr} {}  // Special construct to denote the absence of data
    pointer(const node& subnode) : size_{subnode.size()}, subnode{&subnode} {}
    pointer(const dna& data) : size_{1}, data{&data} {}

    operator std::size_t() const noexcept { return detail::hash(reinterpret_cast<std::size_t>(void_pointer), size_); }

    bool operator==(const pointer& other) const noexcept {
      return size_ == other.size_ && void_pointer == other.void_pointer;
    }

    // Accessor
    auto operator[](std::size_t index) const -> const dna& {
      assert(index < this->size());
      if (this->is_leaf()) return this->get_leaf();
      else return this->get_node()[index];
    }

    auto empty() const noexcept -> bool { return size_ == 0; }
    auto size() const noexcept -> std::size_t { return size_; }
    auto is_leaf() const noexcept -> bool { return size_ == 1; }
    
    auto get_leaf() const -> const dna& { assert(is_leaf() && "Trying to interpret a non-leaf node as a leaf."); return *data; }
    auto get_node() const -> const node& { assert(!is_leaf() && "Trying to interpret a leaf node as a non-leaf."); return *subnode; }

    friend auto operator<<(std::ostream& os, const pointer& p) -> std::ostream& {
      if (p.empty()) return os << "empty";
      else if (p.is_leaf()) return os << "leaf: " << *p.data;
      else return os << "node: " << p.subnode;
    }

  private:
    // TODO: Add annotations for similarity transforms
    std::size_t size_;
    
    union {
      const void* void_pointer;
      const node* subnode;
      const dna* data;
    };
  };

private:
  pointer left, right;
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
  using pointer = typename node::pointer;

  template<typename Iterable>
  shared_tree(Iterable& data);

  auto left_leaf() const -> const dna& { return root.get_node().left_leaf(); }
  auto right_leaf() const -> const dna& { return root.get_node().right_leaf(); }
  auto left_node() const -> const node& { return root.get_node().left_node(); }
  auto right_node() const -> const node& { return root.get_node().right_node(); }

  // Number of elements stored
  auto size() const -> std::size_t { return nodes.size(); }

  // Length of the data sequence
  auto length() const -> std::size_t { return root.size(); }

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
    auto operator!=(const iterator& other) { return &parent != &other.parent || index != other.index; }

  private:
    shared_tree& parent;
    std::size_t index;
  };

  class const_iterator {
  public:
    const_iterator(const shared_tree& parent, std::size_t index = 0) : parent{parent}, index{index} {}
    auto& operator*() { return parent[index]; }
    auto& operator++() { ++index; return *this; }
    auto operator++(int) { auto temp = *this; ++index; return temp; }
    auto& operator--() { --index; return *this; }
    auto operator--(int) { auto temp = *this; --index; return temp; }
    auto operator!=(const const_iterator& other) { return &parent != &other.parent || index != other.index; }

  private:
    const shared_tree& parent;
    std::size_t index;
  };

  auto begin() { return iterator{*this, 0}; }
  auto end() { return iterator{*this, root.size()}; }
  auto cbegin() { return const_iterator{*this, 0}; }
  auto cend() { return const_iterator{*this, root.size()}; }

private:
  pointer root;
  std::unordered_set<node> nodes;
  std::unordered_set<dna> leaves;
};

/**
 *  Constructs a shared binary tree from a vector of data, using spatial
 *  subdivision for common subtree merging.
 */
template<typename Iterable>
shared_tree::shared_tree(Iterable& data) {
    nodes.reserve(data.size()/2 + data.size()%2);
    std::vector<pointer> previous_layer;
    previous_layer.reserve(data.size()/2 + data.size()%2);

    for (const auto& element : data) {
      auto insertion = leaves.emplace(element);
      auto& canonical_leaf = *(insertion.first);
      previous_layer.emplace_back(canonical_leaf);
    }

    while (previous_layer.size() > 1) {
      std::vector<pointer> next_layer;
      next_layer.reserve(previous_layer.size()/2 + previous_layer.size()%2);

      for (const auto& [left, right] : pairwise(previous_layer)) {
        auto created_node = node{left, right};
        auto& canonical_node = *(nodes.emplace(created_node).first);
        next_layer.emplace_back(canonical_node);
      }
      if (previous_layer.size() % 2) {
        auto created_node = node{previous_layer.back()};
        auto& canonical_node = *(nodes.emplace(created_node).first);
        next_layer.push_back(canonical_node);
      }

      previous_layer = std::move(next_layer);
    }

    root = pointer{previous_layer.front()};
  }