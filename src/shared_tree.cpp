/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 */
 
#include <algorithm>
#include <fstream>

 #include "shared_tree.h"

/******************************************************************************
 *  Annotated pointer to a node in the tree.
 */
auto pointer::leaf() const -> dna {
  assert(!empty() && "Trying to access null pointer.");
  assert(is_leaf() && "Trying to interpret a non-leaf node as a leaf.");
  return dna{information};
}

auto pointer::index() const -> std::uint64_t {
  assert(!empty() && "Trying to access null pointer.");
  assert(!is_leaf() && "Trying to interpret a leaf node as a non-leaf.");
  return information;
}

auto pointer::to_ullong() const noexcept -> unsigned long long {
  return information
  | static_cast<unsigned long long>(mirror_bit) << 60
  | static_cast<unsigned long long>(transpose_bit) << 61
  | static_cast<unsigned long long>(leaf_bit) << 62;
}


/******************************************************************************
 *  Node in a binary tree.
 */
/**
 * Two nodes are considered equal if they can be transformed into one another
 * through only mirroring and/or transposition.
 */
bool node::operator==(const node& other) const noexcept {
  const auto transpose = std::array{children[0].transposed(), children[1].transposed()};
  const auto mirror = std::array{children[1].mirrored(), children[0].mirrored()};
  const auto both = std::array{children[1].mirrored().transposed(), children[0].mirrored().transposed()};

  const bool match_exact = children[0] == other.children[0] && children[1] == other.children[1];
  const bool match_transposed = transpose[0] == other.children[0] && transpose[1] == other.children[1];
  const bool match_mirrored = mirror[0] == other.children[0] && mirror[1] == other.children[1];
  const bool match_both = both[0] == other.children[0] && both[1] == other.children[1];

  return match_exact || match_transposed || match_mirrored || match_both;
}

/**
 * Determines which transformations need to be applied to obtain the canonical
 * node. Returns them in a pair, first member referring to mirroring, second to
 * transposition.
 * Precondition: other must be similar to *this.
 */
auto node::transformations(const node& other) const noexcept -> std::pair<bool, bool> {
  const auto transpose = std::array{children[0].transposed(), children[1].transposed()};
  const auto mirror = std::array{children[1].mirrored(), children[0].mirrored()};

  const bool match_exact = children[0] == other.children[0] && children[1] == other.children[1];
  const bool match_transposed = transpose[0] == other.children[0] && transpose[1] == other.children[1];
  const bool match_mirrored = mirror[0] == other.children[0] && mirror[1] == other.children[1];

  if (match_exact) return {false, false};
  if (match_transposed) return {false, true};
  if (match_mirrored) return {true, false};
  else return {true, true};
}


/******************************************************************************
 *  Genome, constructed as a balanced shared tree.
 */
/**
 *  Returns the number of child leaves the node pointed to has.
 */
auto shared_tree::children(pointer parent) const -> std::size_t {
  if (parent.is_leaf()) return 1;
  if (parent.empty()) return 0;
  auto node = access(parent);
  return children(node.left()) + children(node.right());
}

/**
 * Constructs a histogram of the number of references to each node, and stores
 * it in .csv format at the provided location.
 */
void shared_tree::histogram(std::filesystem::path path) const {
  auto frequencies = std::vector<unsigned long long>(nodes.size(), 0);
  auto file = std::ofstream{path};

  for (const auto& node : nodes) {
    if (!node.left().is_leaf() && !node.left().empty()) ++frequencies[node.left().index()-1];
    if (!node.right().is_leaf() && !node.right().empty()) ++frequencies[node.right().index()-1];
  }

  std::sort(frequencies.begin(), frequencies.end(), std::greater<>());

  for (const auto& frequency : frequencies)
    file << frequency << ',';
}

/**
 *  Constructs an iterator over the shared tree.
 *  A nullptr argument for the root indicates an end iterator.
 */
shared_tree::iterator::iterator(const std::vector<node>& nodes, pointer root) : nodes{nodes} {
  if (root != nullptr) {
    stack.emplace_back(root);
    transposed.emplace_back(root.is_transposed());
    next_leaf();
  }
}

/**
 * Returns the DNA strand stored in the current leaf.
 * Mirrors and transposes it, if necessary.
 */
auto shared_tree::iterator::operator*() noexcept -> dna {
  auto result = stack.back().leaf();
  if (mirrored) result = result.mirrored();
  if (transposed.back()) result = result.transposed();
  return result;
}

/**
 *  Advances the iterator to the first node on the right of the current node.
 *  To achieve this, pops the last-found leaf off the stack and recurses on the
 *  remaining stack members.
 */
auto shared_tree::iterator::operator++() -> iterator& {
  stack.pop_back();
  transposed.pop_back();
  next_leaf();
  return *this;
}

/**
 *  Iterates over the tree nodes until a leaf is found, or the stack is empty.
 *  For each encountered node, determines if it is a leaf. If not, its children
 *  are added to the stack and we apply the same procedure to them, starting
 *  with the left child.
 */
void shared_tree::iterator::next_leaf() {
  while (!stack.empty()) {
    auto top = stack.back();
    mirrored = top.is_mirrored();

    if (top.is_leaf()) return;

    auto node = access(top);
    stack.pop_back();
    transposed.pop_back();

    if (!mirrored) {
      if (auto right = node.right(); !right.empty()) {
        stack.emplace_back(right);
        transposed.emplace_back(right.is_transposed() ^ top.is_transposed());
      }
      if (auto left = node.left(); !left.empty()) {
        stack.emplace_back(left);
        transposed.emplace_back(left.is_transposed() ^ top.is_transposed());
      }
    } else {
      if (auto left = node.left(); !left.empty()) {
        stack.emplace_back(left);
        transposed.emplace_back(left.is_transposed() ^ top.is_transposed());
      }
      if (auto right = node.right(); !right.empty()) {
        stack.emplace_back(right);
        transposed.emplace_back(right.is_transposed() ^ top.is_transposed());
      }
    }
  }
}

/**
 *  Access function identical to the one in class shared_tree itself.
 */
auto shared_tree::iterator::access(pointer pointer) const -> const node& {
  const auto index = pointer.index();
  return nodes[index-1];
}

/**
 *  Since the pointers stored inside nodes are not sufficient without a
 *  reference to the data structure holding them, the accessor logic must
 *  reside in the shared tree.
 */
auto shared_tree::operator[](std::size_t index) const -> dna {
  auto current = root;
  bool mirror = false;
  bool transpose = false;

  while (!current.is_leaf()) {
    mirror = current.is_mirrored();
    transpose ^= current.is_transposed();
    const auto& node = access(current);

    if (!mirror) {
      const auto left_size = children(node.left());
      if (index < left_size) {
        current = node.left();
      } else {
        index -= left_size;
        current = node.right();
      }
    } else {
      const auto right_size = children(node.right());
      if (index < right_size) {
        current = node.right();
      } else {
        index -= right_size;
        current = node.left();
      }
    }
  }

  auto leaf = current.leaf();
  transpose ^= current.is_transposed();
  mirror = current.is_mirrored();
  if (transpose) leaf = leaf.transposed();
  if (mirror) leaf = leaf.mirrored();
  return leaf;
}

/**
 *  Returns the node pointed to by the passed pointer.
 *  Accessing a nullptr results in undefined behaviour.
 */
auto shared_tree::access(pointer pointer) const -> const node& {
  const auto index = pointer.index();
  return nodes[index-1];
}