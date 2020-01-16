/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 *  Specialised for balanced trees, to improve compression speed and ratio for
 *  those trees. This does mean that unbalanced trees are not supported.
 */

#include <limits>

#include "balanced_shared_tree.h"

#include "fasta_reader.h"

/****************************************************************************
 * class node:
 *  Node type representing inner nodes in the tree.
 *  Consists of two pointers to nodes or leaves one level down in the tree.
 */
/**
 * Two nodes compare equal when they can be transformed into one another
 * through any possible combination of similarity transforms.
 */
bool detail::node::operator==(const node& other) const noexcept {
  return children == other.children
  || mirrored().children == other.children
  || transposed().children == other.children
  || inverted().children == other.children;
}

/**
 * Determines the transformations required to obtain other from the current
 * node. First boolean corresponds to mirroring, second to transposition.
 * Precondition: other must be similar to this.
 */
auto detail::node::transformations(const node& other) const noexcept -> std::pair<bool, bool> {
  if (children == other.children) return {false, false};
  if (mirrored().children == other.children) return {true, false};
  if (transposed().children == other.children) return {false, true};
  else return {true, true};
}


/******************************************************************************
 * class balanced_shared_tree:
 *  Shared binary tree class that exploits structural properties of balanced
 *  trees to store its directed acyclic graph representation more effectively.
 */
/**
 * Constructs a balanced_shared_tree from a FASTA formatted file.
 */
balanced_shared_tree::balanced_shared_tree(std::filesystem::path path) {
  constexpr auto segment_size = (1u<<29);
  auto data = fasta_reader{path};
  auto constructor = tree_constructor{*this};
  auto subroots = std::vector<pointer>{};

  for (auto segment : chunks(data, segment_size))
    constructor.reduce(segment);
  constructor.reduce();
  root = constructor.root();
}

/**
 * Returns the total number of nodes in the tree.
 * Leafs are excluded.
 */
auto balanced_shared_tree::node_count() const -> std::size_t {
  auto sum = 0ull;
  for (const auto& layer : nodes) sum += layer.size();
  return sum;
}

/**
 * Accesses the leaf pointed to by <pointer>.
 * Applies the correct transforms to the canonical leaf to obtain the result.
 */
auto balanced_shared_tree::access_leaf(pointer pointer) const -> dna {
  auto result = leaves[pointer.index()];
  if (pointer.is_mirrored()) result = result.mirrored();
  if (pointer.is_transposed()) result = result.transposed();
  return result;
}

/**
 * Accesses the node in layer <layer> pointed to by <pointer>.
 * Returned by value since since the nodes must be immutable anyway.
 */
auto balanced_shared_tree::access_node(std::size_t layer, pointer pointer) const -> node {
  return nodes[layer][pointer.index()];
}

/**
 * Returns the number of children contained in the subtree referenced by
 * <pointer> in <layer>.
 */
auto balanced_shared_tree::children(std::size_t layer, pointer pointer) const -> std::size_t {
  if (pointer == nullptr) return 0;
  const auto node = access_node(layer, pointer);
  const auto left = node.left();
  const auto right = node.right();
  if (layer == 0) return !left.empty() + !right.empty();
  else return children(layer-1, left) + children(layer-1, right);
}

auto balanced_shared_tree::operator[](std::uint64_t index) const -> dna {
  auto current = root;
  auto layer = nodes.size() - 1;

  auto combine = [](auto next, auto current) {
    const auto mirror = next.is_mirrored() != current.is_mirrored();
    const auto transpose = next.is_transposed() != current.is_transposed();
    return pointer{next.index(), mirror, transpose};
  };

  auto descend = [&](auto left, auto right) {
    const auto size = children(--layer, left);
    if (index < size) return combine(left, current);
    index -= size;
    return combine(right, current);
  };

  while (layer > 0) {
    const auto& node = access_node(layer, current);
    const auto left = node.left();
    const auto right = node.right();

    if (current.is_mirrored()) current = descend(right, left);
    else current = descend(left, right);
  }

  const auto& node = access_node(0, current);
  const auto left = node.left();
  const auto right = node.right();

  auto descend_leaf = [&](auto left, auto right) {
    const auto size = (bool)left;
    if (index < size) return combine(left, current);
    index -= size;
    return combine(right, current);
  };

  if (current.is_mirrored()) current = descend_leaf(right, left);
  else current = descend_leaf(left, right);
  return access_leaf(current);
}

/**
 * Adds a node to the specified layer.
 * Precondition: no similar nodes are already present in this layer. 
 */
void balanced_shared_tree::emplace_node(std::size_t layer, node node) {
  nodes[layer].emplace_back(node);
}

/**
 * Adds a leaf to the leaf layer.
 * Precondition: the leaf is in its canonical shape and not yet present.
 */
void balanced_shared_tree::emplace_leaf(dna leaf) {
  leaves.emplace_back(leaf);
}


/******************************************************************************
 * class balanced_shared_tree::iterator:
 *  Iterator over the tree.
 */
balanced_shared_tree::iterator::iterator(balanced_shared_tree& parent, std::size_t layer, detail::pointer root) : parent{parent} {
  if (root) {
    stack.emplace_back(layer, root);
    next_leaf();
  }
}

/**
 * Returns the DNA strand stored in the current leaf.
 * Mirrors and transposes it, if necessary.
 * Precondition: current stack top points to a leaf.
 */
auto balanced_shared_tree::iterator::operator*() const noexcept -> dna {
  auto top = stack.back().current;
  return parent.access_leaf(top);
}

/**
 *  Advances the iterator to the first leaf on the right of the current leaf.
 *  To achieve this, pops the last-found leaf off the stack and recurses on the
 *  remaining stack members until a leaf is found.
 */
auto balanced_shared_tree::iterator::operator++() -> iterator& {
  stack.pop_back();
  next_leaf();
  return *this;
}

/**
 *  Iterates over the tree nodes until a leaf is found, or the stack is empty.
 *  For each encountered node, determines if it is a leaf. If not, its children
 *  are added to the stack and we apply the same procedure to them, starting
 *  with the left child.
 */
void balanced_shared_tree::iterator::next_leaf() {
  while (!stack.empty()) {
    auto status = stack.back();
    // Leaf level is indicated as the maximum std::size_t value
    if (status.layer == std::numeric_limits<std::size_t>::max()) return;

    auto top = status.current;
    auto node = parent.access_node(status.layer, status.current);
    stack.pop_back();

    // Apply mirroring and transposition, if necessary, and save the resulting
    // state on the stack.
    auto stack_push = [&](auto next) {
        auto mirror = next.is_mirrored() != top.is_mirrored();
        auto transpose = next.is_transposed() != top.is_transposed();
        auto updated_pointer = pointer{next.index(), mirror, transpose};
        stack.emplace_back(status.layer - 1, updated_pointer);
    };

    if (top.is_mirrored()) {
      if (node.left()) stack_push(node.left());
      if (node.right()) stack_push(node.right());
    } else {
      if (node.right()) stack_push(node.right());
      if (node.left()) stack_push(node.left());
    }
  }
}

/******************************************************************************
 * class tree_constructor:
 *  Helper class in construction of a balanced shared tree.
 *  Contains the maps used to link nodes to pointers or leaves.
 */
tree_constructor::tree_constructor(balanced_shared_tree& parent)
: parent{parent} {}

/**
 * Constructs and emplaces a node inside the tree during its construction.
 * A pointer to the node is stored in <next_layer>.
 */
void tree_constructor::emplace_node(pointer left, pointer right) {
  auto created_node = node{left, right};
  auto insertion = nodes[layer].emplace(created_node, parent.node_count(layer));
  auto canonical_node = (*insertion.first).first;
  auto index = (*insertion.first).second;

  if (insertion.second) {
    parent.emplace_node(layer, created_node);
    next_layer.emplace_back(index, false, false);
  } else {
    auto transform = created_node.transformations(canonical_node);
    next_layer.emplace_back(index, transform.first, transform.second);
  }
}

/**
 * Constructs and emplaces a leaf inside the tree during its construction.
 * A pointer to the leaf is stored in <current_layer>, in preparation for
 * further reduction.
 */
void tree_constructor::emplace_leaf(const dna& leaf) {
  auto [canonical_leaf, mirror, transpose] = leaf.canonical();
  auto insertion = leaves.emplace(canonical_leaf, parent.leaf_count());
  auto index = (*insertion.first).second;

  if (insertion.second) parent.emplace_leaf(canonical_leaf);
  current_layer.emplace_back(index, mirror, transpose);
}

/**
 * Reduces the current layer by constructing nodes out of the pointers it
 * contains, and emplacing those nodes in the correct layers and maps, if
 * necessary.
 */
void tree_constructor::reduce_once() {
  next_layer.reserve(current_layer.size()/2 + current_layer.size()%2);
  if (parent.depth() <= layer) {
    parent.add_layer();
    nodes.emplace_back();
  }

  foreach_pair(current_layer,
    [&](const auto& left, const auto& right) { emplace_node(left, right); },
    [&](const auto& last) { emplace_node(last); }
  );

  current_layer = std::move(next_layer);
  next_layer.clear();
  ++layer;
}


/**
 * Combines the constructed subtrees into a single tree by further reduction.
 * Precondition: <layer> is valued at the correct index for the subroots.
 */
void tree_constructor::reduce() {
  current_layer = roots;
  roots.clear();
  while (current_layer.size() > 1) {
    reduce_once();
  } 
  roots.emplace_back(current_layer.front());
}