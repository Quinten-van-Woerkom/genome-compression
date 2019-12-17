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

#include "utility.h"

/**
 *  Node in a non-owning tree. Children are either other nodes or leaf nodes,
 *  containing the actual data stored.
 */
class node {
public:
  class pointer;

  node(pointer left, pointer right) : left_{left}, right_{right} {}
  node(pointer left) : left_{left}, right_{nullptr} {}
  node() : left_{nullptr}, right_{nullptr} {}

  auto left_leaf() const -> const dna& { assert(left_.is_leaf()); return left_.get_leaf(); }
  auto right_leaf() const -> const dna& { assert(right_.is_leaf()); return right_.get_leaf(); }
  auto left_node() const -> const node& { assert(!left_.is_leaf()); return left_.get_node(); }
  auto right_node() const -> const node& { assert(!right_.is_leaf()); return right_.get_node(); }

  auto left() const -> pointer { return left_; }
  auto right() const -> pointer { return right_; }
  
  auto hash() const noexcept -> std::size_t { return detail::hash(left_, right_); }
  auto size() const noexcept -> std::size_t { return left_.size() + right_.size(); }
  auto full() const noexcept -> bool { return !left_.empty() && !right_.empty(); }
  auto empty() const noexcept -> bool { return left_.empty() && right_.empty(); }
  void clear() noexcept { left_ = nullptr; right_ = nullptr; }

  template<typename... Args>
  auto left(Args&&... args) {
    left_ = pointer{args...};
    return left_;
  }

  template<typename... Args>
  auto right(Args&&... args) {
    right_ = pointer{args...};
    return right_;
  }

  // Tries to add a new pointer in the first empty position, if there is any.
  // Returns false if both left and right positions are filled.
  template<typename... Args>
  auto emplace(Args&&... args) noexcept -> bool {
    if (left_.empty()) {
      left(args...);
      return true;
    } else if (right_.empty()) {
      right(args...);
      return true;
    } else return false;
  }

  // For now, only mirror symmetry exists.
  // In the future, more similarities might be added.
  auto similarities(const node& other) const {
    assert(*this == other);
    if (left_ == other.right_ && right_ == other.left_) return true;
    else return false;
  }

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

    auto empty() const noexcept -> bool { return size_ == 0 && void_pointer == nullptr; }
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
  using pointer = typename node::pointer;

  shared_tree(const shared_tree&) = delete; // Ensures that factory functions are used for in-place construction
  shared_tree() : root{nullptr} {};

  shared_tree(shared_tree&& other)
    : root{std::move(other.root)}, nodes{std::move(other.nodes)},
      leaves{std::move(other.leaves)}
  {}

  template<typename Iterable>
  static auto create_balanced_layerwise(Iterable&& data) -> shared_tree;

  template<typename Iterable>
  static auto create_balanced_pairwise(Iterable&& data) -> shared_tree;

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
 *  Reduces the data into nodes layer-by-layer; this requires repeated linear
 *  traversal of data, which is fast but requires a lot of memory.
 */
template<typename Iterable>
auto shared_tree::create_balanced_layerwise(Iterable&& data) -> shared_tree {
  auto result = shared_tree{};
  result.nodes.reserve(data.size()/2 + data.size()%2);
  std::vector<pointer> previous_layer;
  previous_layer.reserve(data.size()/2 + data.size()%2);

  for (const auto& element : data) {
    auto insertion = result.leaves.emplace(element);
    auto& canonical_leaf = *(insertion.first);
    previous_layer.emplace_back(canonical_leaf);
  }
  
  while (previous_layer.size() > 1) {
    std::vector<pointer> next_layer;
    next_layer.reserve(previous_layer.size()/2 + previous_layer.size()%2);

    for (const auto& [left, right] : pairwise(previous_layer)) {
      auto created_node = node{left, right};
      auto insertion = result.nodes.emplace(created_node);
      auto& canonical_node = *(insertion.first);
      next_layer.emplace_back(canonical_node);
    }
    if (previous_layer.size() % 2) {
      auto created_node = node{previous_layer.back()};
      auto insertion = result.nodes.emplace(created_node);
      auto& canonical_node = *(insertion.first);
      next_layer.emplace_back(canonical_node);
    }

    previous_layer = std::move(next_layer);
  }

  result.root = pointer{previous_layer.front()};
  return result;
}

/**
 *  Constructs a shared binary tree from a vector of data, using spatial
 *  subdivision for common subtree merging.
 *  Reduces the data into nodes pair-by-pair; this reduces memory requirements
 *  and allows for a single traversal of the original dataset to be sufficient.
 */
// template<typename Iterable>
// auto shared_tree::create_balanced_pairwise(Iterable&& data) -> shared_tree {
//   auto result = shared_tree();
//   auto leaves = std::vector<const dna*>{};
//   auto pointers = std::vector<std::vector<pointer>>{};
//   pointers.emplace_back();

//   for (const auto& element : data) {
//     auto insertion = result.leaves.emplace(element);
//     const auto& canonical_leaf = *(insertion.first);
//     leaves.emplace_back(&canonical_leaf);
    
//     if (leaves.size() == 2) {
//       auto created_node = node{*leaves[0], *leaves[1]};
//       auto leaf_insertion = result.nodes.emplace(created_node);
//       auto& canonical_node = *(leaf_insertion.first);
//       pointers[0].emplace_back(pointer{canonical_node});

//       for (auto level = 0u; pointers[level].size() == 2 && level < pointers.size(); ++level) {
//         auto new_node = node{pointers[level][0], pointers[level][1]};
//         auto node_insertion = result.nodes.emplace(new_node);
//         auto& canonical_new_node = *(node_insertion.first);
        
//         if (level+1 == pointers.size()) {
//           pointers.emplace_back();
//         }
//         pointers[level+1].emplace_back(pointer{canonical_new_node});
//         pointers[level].clear();
//       }

//       leaves.clear();
//     }
//   }
  
//   if (leaves.size() == 1) {
//     std::cout << "Cleaning up singular leave\n";
//     auto created_node = node{*leaves[0]};
//     auto insertion = result.nodes.emplace(created_node);
//     auto& canonical_node = *(insertion.first);
//     pointers[0].emplace_back(pointer{canonical_node});
//   }

//   for (auto level = 0u; level < pointers.size(); ++level) {
//     std::cout << "Cleaning up singular nodes\n";
//     if (pointers[level].size() == 1) {
//       auto new_node = node{pointers[level][0]};
//       auto node_insertion = result.nodes.emplace(new_node);
//       auto& canonical_new_node = *(node_insertion.first);
      
//       if (level+1 == pointers.size() && pointers[level].size() == 2) pointers.emplace_back();
//       pointers[level+1].emplace_back(pointer{canonical_new_node});
//       pointers[level].clear();
//     } else {
//       auto new_node = node{pointers[level][0], pointers[level][1]};
//       auto node_insertion = result.nodes.emplace(new_node);
//       auto& canonical_new_node = *(node_insertion.first);
      
//       if (level+1 == pointers.size()) pointers.emplace_back();
//       pointers[level+1].emplace_back(pointer{canonical_new_node});
//       pointers[level].clear();
//     }
//   }

//   std::cout << "Creating root\n";
//   result.root = [&]() {
//     if (pointers.back().size() == 1) return pointers.back()[0];
//     else {
//       auto& canonical_root = *(result.nodes.emplace(pointers.back()[0], pointers.back()[1]).first);
//       return pointer{canonical_root};
//   }}();

//   return result;
// }

template<typename Iterable>
auto shared_tree::create_balanced_pairwise(Iterable&& data) -> shared_tree {
  auto result = shared_tree();
  auto nodes = std::vector<node>{};
  nodes.emplace_back();

  for (const auto& element : data) {
    auto insertion = result.leaves.emplace(element);
    const auto& leaf = *insertion.first;
    nodes[0].emplace(leaf);
    
    for (auto level = 0u; level < nodes.size(); ++level) {
      if (nodes[level].full()) {
        auto insertion = result.nodes.emplace(nodes[level]);
        const auto& node = *insertion.first;
        if (level+1 == nodes.size()) nodes.emplace_back(pointer{node});
        else nodes[level+1].emplace(pointer{node});
        nodes[level].clear();
      }
    }
  }

  for (auto level = 0u; level < nodes.size()-1; ++level) {
    auto insertion = result.nodes.emplace(nodes[level]);
    const auto& node = *insertion.first;
    nodes[level+1].emplace(pointer{node});
  }

  result.root = nodes.back();
  return result;
}