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
balanced_shared_tree::balanced_shared_tree(fasta_reader file) {
  auto constructor = tree_constructor{*this};
  root = constructor.reduce(file);
}

balanced_shared_tree::balanced_shared_tree(std::vector<dna>& data) {
  auto constructor = tree_constructor{*this};
  root = constructor.reduce(data);
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
 * Accesses the first or second leaf located in the node pointed to by
 * <pointer>, depending on <index>.
 */
auto balanced_shared_tree::access_leaf(std::size_t index, pointer current) const -> dna {
  auto parent_node = access_node(0, current);
  const auto left = parent_node.left();
  const auto right = parent_node.right();

  auto combine = [](auto next, auto current) {
    const auto mirror = next.is_mirrored() != current.is_mirrored();
    const auto transpose = next.is_transposed() != current.is_transposed();
    return pointer{next.index(), mirror, transpose};
  };

  auto descend_leaf = [&](auto left, auto right) {
    const auto size = (bool)left;
    if (index < size) return combine(left, current);
    index -= size;
    return combine(right, current);
  };
  
  if (current.is_mirrored()) current = descend_leaf(right, left);
  else current = descend_leaf(left, right);
  return current.leaf();
}

/**
 * Accesses the node in layer <layer> pointed to by <pointer>.
 * Returned by value since since the nodes must be immutable anyway.
 */
auto balanced_shared_tree::access_node(std::size_t layer, pointer pointer) const -> node {
  assert(layer < nodes.size());
  if (pointer.index() >= nodes[layer].size())
    std::cerr << "Pointer index out of bounds: " << pointer.index() << " for a size of " << nodes[layer].size() << " in layer " << layer << '\n';
  return nodes[layer][pointer.index()];
}

/**
 * Returns the number of children contained in the subtree referenced by
 * <pointer> in <layer>. To determine this, we just traverse the tree without
 * considering mirroring or transposition, as those do not alter the number of
 * children a node has.
 */
auto balanced_shared_tree::children(std::size_t layer, pointer pointer) const -> std::size_t {
  if (pointer.empty()) return 0;
  const auto node = access_node(layer, pointer);
  const auto left = node.left();
  const auto right = node.right();
  if (layer == 0) return !left.empty() + !right.empty();
  else return children(layer-1, left) + children(layer-1, right); 
}

auto balanced_shared_tree::operator[](std::uint64_t index) const -> dna {
  auto current = root;

  auto combine = [](auto next, auto current) {
    const auto mirror = next.is_mirrored() != current.is_mirrored();
    const auto transpose = next.is_transposed() != current.is_transposed();
    return pointer{next.index(), mirror, transpose};
  };

  auto descend = [&](auto layer, auto left, auto right) {
    const auto size = children(layer, left);
    if (index < size) return combine(left, current);
    index -= size;
    return combine(right, current);
  };

  for (auto layer = nodes.size()-1; layer > 0; --layer) {
    const auto& node = access_node(layer, current);
    const auto left = node.left();
    const auto right = node.right();

    if (current.is_mirrored()) current = descend(layer-1, right, left);
    else current = descend(layer-1, left, right);
  }

  return access_leaf(index, current);
}

/**
 * Adds a node to the specified layer.
 * Precondition: no similar nodes are already present in this layer. 
 */
void balanced_shared_tree::emplace_node(std::size_t layer, node node) {
  nodes[layer].emplace_back(node);
}

/**
 * Creates and stores a histogram for each layer in the tree, storing them as
 * lines in a .csv file.
 */
void balanced_shared_tree::histogram(std::filesystem::path path) const {
  auto file = std::ofstream{path};
  for (auto layer = 1u; layer < nodes.size(); ++layer) {
    auto frequencies = std::vector<unsigned long long>(nodes[layer-1].size(), 0);

    // TODO: this is bugged, fix it
    for (const auto& node : nodes[layer]) {
      // if (!node.left().empty() && node.left().index() < frequencies.size())
      if (!node.left().empty())
        ++frequencies[node.left().index()];
      // if (!node.right().empty() && node.right().index() < frequencies.size())
      if (!node.right().empty())
        ++frequencies[node.right().index()];
    }

    std::sort(frequencies.begin(), frequencies.end(), std::greater<>());

    for (const auto& frequency : frequencies)
      file << frequency << ',';
    file << '\n';
  }
}

/**
 * Prints all layers to the output stream, in order from the root to the
 * leaves.
 */
void balanced_shared_tree::print_unique(std::ostream& os) const {
  os << "Root (" << nodes.back().size() << "): ";
  for (auto node : nodes.back()) os << node << ", ";
  os << '\n';

  for (int i = nodes.size()-2; i > 0; --i) {
    os << "Layer " << i << " (" << nodes[i].size() << "): ";
    for (auto node : nodes[i]) {
      if (node.left() >= nodes[i-1].size() || node.right() >= nodes[i-1].size())
        os << node << ", ";
    }
    os << '\n';
  }

  // os << "Leaves (" << nodes.front().size() << "): ";
  // for (auto node : nodes.front()) os << node.left().leaf() << ' ' << node.right().leaf() << ", ";
  // os << '\n';
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
  return top.leaf();
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
void tree_constructor::emplace(
  std::vector<pointer>& layer, std::size_t layer_index,
  pointer left, pointer right)
{
  auto created_node = node{left, right};
  auto insertion = nodes[layer_index].emplace(created_node, parent.node_count(layer_index));
  auto canonical_node = (*insertion.first).first;
  auto index = (*insertion.first).second;

  if (insertion.second) {
    parent.emplace_node(layer_index, created_node);
    layer.emplace_back(index, false, false);
  } else {
    auto transform = created_node.transformations(canonical_node);
    layer.emplace_back(index, transform.first, transform.second);
  }
}