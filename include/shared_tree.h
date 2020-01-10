/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 */

#pragma once

#include <cassert>
#include <optional>
#include <stack>
#include <unordered_map>
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
  // And index of 0 is used to represent a nullptr, and absence of data
  pointer(std::size_t index) : leaf_flag{false}, data{index} {}
  pointer(dna dna) : leaf_flag{true}, data{dna.to_ullong()} {}
  pointer(std::nullptr_t) : leaf_flag{false}, data{0} {}

  operator std::size_t() const noexcept { return detail::hash(leaf_flag, data); }  // TODO: Check if this is really necessary.

  bool operator==(const pointer& other) const noexcept {
    return leaf_flag == other.leaf_flag && data == other.data;
  }

  auto is_leaf() const noexcept -> bool { return leaf_flag; }
  auto empty() const noexcept -> bool { return !is_leaf() && data == 0; }
  auto leaf() const -> dna;
  auto index() const -> std::size_t;

private:
  bool leaf_flag : 1;
  std::size_t data : 60;
};

static auto operator<<(std::ostream& os, const pointer& p) -> std::ostream& {
  if (p.empty()) return os << "empty";
  else if (p.is_leaf()) return os << "leaf: " << p.leaf();
  else return os << "node: " << p.index();
}

/******************************************************************************
 *  Node in a shared tree. Children are either other nodes or leaf nodes,
 *  containing the actual data stored.
 */
class node {
public:
  node(pointer left, pointer right) : left_{left}, right_{right} {}
  node(pointer left) : left_{left}, right_{nullptr} {}
  
  auto left() const -> const pointer& { return left_; }
  auto right() const -> const pointer& { return right_; }
  
  auto hash() const noexcept -> std::size_t { return detail::hash(left_, right_); }
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
  auto width() const -> std::size_t { return children(root); }
  auto children(pointer parent) const -> std::size_t;

  auto operator[](std::size_t index) const -> dna;
  
  void print_unique() const {
    for (const auto& node : nodes)
        std::cout << node << '\n';
  }

  struct iterator {
    iterator(const std::vector<node>& nodes, pointer root) : nodes{&nodes}, current{root} {
      if (root != pointer{nullptr}) {
        stack.emplace_back(root);
        ++(*this);
      }
    }

    auto operator*() const { return current.leaf(); }
    auto operator++() -> iterator&;
    auto operator!=(const iterator& other) const { return nodes != other.nodes || current != other.current; }

    auto access(pointer pointer) const -> const node&;

    const std::vector<node>* nodes;
    pointer current;
    std::vector<pointer> stack;
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

private:
  auto access(pointer pointer) const -> const node&;

  template<typename... Args>
  void emplace_node(std::vector<pointer>& layer, Args&&... args);

  template<typename Iterable>
  auto reduce_once(Iterable&& layer) -> std::vector<pointer>;

  template<typename Iterable>
  auto reduce(Iterable&& layer) -> pointer;

  pointer root;
  std::unordered_map<node, std::size_t> indices;
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
  auto created_node = node{pointer{args}...};
  auto insertion = indices.try_emplace(created_node, nodes.size()+1);
  if (insertion.second) nodes.emplace_back(created_node);
  auto index = (*insertion.first).second;
  layer.emplace_back(index);
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