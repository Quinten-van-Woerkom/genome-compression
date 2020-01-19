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
 * class pointer:
 *  Pointer type representing references to another node or to a leaf node.
 */
/**
 * Helper function that returns the width of the address space corresponding to
 * the given segment.
 */
constexpr auto address_space(std::size_t segment) noexcept {
  return 1ull << detail::pointer::address_bits[segment];
}

/**
 * Constructor from another pointer. Transformations can be applied to the
 * pointer, if necessary.
 * Note that a transformed nullptr must also be a nullptr, so that bits must
 * be set in this special case.
 */
detail::pointer::pointer(const pointer& other, bool mirror, bool transpose)
: data{other.data},
  mirror{mirror != other.mirror || other == nullptr},
  transpose{transpose != other.transpose || other == nullptr},
  segment{other.segment} {}

/**
 * Constructor from an index with given similarity transforms.
 * Stores the index in (segment, offset) format, so that smaller indices can
 * be represented as shorter pointers.
 */
detail::pointer::pointer(std::size_t index, bool mirror, bool transpose)
: data{index}, mirror{mirror}, transpose{transpose}, segment{0} {
  while (data >= address_space(segment)) {
    ++segment;
    data -= address_space(segment-1);
  }
}

/**
 * Constructs a direct pointer referring to a DNA strand. Since the DNA strands
 * are smaller in memory than a pointer, they are stored directly.
 * To keep the pointer-like nature of the data, we store the canonical version
 * of a leaf node, which is the version among the set of similar leaves that
 * has the smallest bit representation.
 */
detail::pointer::pointer(dna leaf, bool add_mirror, bool add_transpose)
: segment{0b11} {
  auto [canonical, mirror_bit, transpose_bit] = leaf.canonical();
  data = canonical.to_ullong();
  mirror = mirror_bit != add_mirror;
  transpose = transpose_bit != add_transpose;
}

/**
 * Constructs a null pointer, indicating an empty subtree.
 * Null pointers have all bits set and are of maximum pointer size. Due to
 * their infrequent occurrence, this bigger size is not an issue, while it
 * significantly simplifies the indexing code.
 * Leaf pointers can be null, too; a fully set bit sequence is not a valid
 * canonical leaf, so there are no collisions with valid DNA strands.
 */
detail::pointer::pointer(std::nullptr_t)
: mirror{true}, transpose{true}, segment{0b11} {
  data ^= ~data;
}


/**
 * Interprets the data as an index pointing to an inner node.
 */
auto detail::pointer::index() const noexcept -> std::size_t {
  assert(!empty());
  auto offset = data;
  if (segment >= 0b11) offset += address_space(0b10);
  if (segment >= 0b10) offset += address_space(0b01);
  if (segment >= 0b01) offset += address_space(0b00);
  return offset;
}

/**
 * Interprets the data contained in the pointer as a canonical leaf with
 * applied transformations.
 */
auto detail::pointer::leaf() const noexcept -> dna {
  assert(!empty());
  auto result = dna{data};
  if (mirror) result = result.mirrored();
  if (transpose) result = result.transposed();
  return result;
}

/**
 * Serializes the pointer to the given output stream in compressed format.
 * Starts with the 4 header bits and the most significant 4 bits, then
 * repeatedly stores the next most significant byte, until no more bytes
 * remain.
 */
void detail::pointer::serialize(std::ostream& os) const {
  int index = address_bits[segment]-4;
  std::uint8_t store = ((data >> index) & 0xf) | mirror << 4 | transpose << 5 | segment << 6;
  os << store;
  for (index -= 8; index >= 0; index -= 8) {
    store = static_cast<std::uint8_t>(data >> index);
    os << store;
  }
}

/**
 * Loads a pointer from an input stream.
 */
auto detail::pointer::deserialize(std::istream& is) -> pointer {
  auto result = pointer{};
  std::uint8_t loaded;
  is.get(reinterpret_cast<char&>(loaded));
  auto index = address_bits[result.segment]-4;
  
  result.segment = (loaded >> 6) & 3u;
  result.transpose = (loaded >> 5) & 1u;
  result.mirror = (loaded >> 4) & 1u;
  result.data = ((std::uint64_t)loaded & 0xf) << index;

  for (index -= 8; index >= 0; index -= 8) {
    is.get(reinterpret_cast<char&>(loaded));
    result.data |= ((std::uint64_t)loaded << index);
  }
  return result;
}


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

/**
 * Serializes the node to the given output stream.
 * First, the left pointer is stored, and then the right pointer.
 */
void detail::node::serialize(std::ostream& os) const {
  children[0].serialize(os);
  children[1].serialize(os);
}

/**
 * Deserializes a node from an input stream.
 * The node is assumed to be stored in compressed format.
 */
auto detail::node::deserialize(std::istream& is) -> node {
  auto left = pointer::deserialize(is);
  auto right = pointer::deserialize(is);
  return node{left, right};
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

  auto descend_leaf = [&](auto left, auto right) {
    const auto size = (bool)left;
    const auto mirror = current.is_mirrored();
    const auto transpose = current.is_transposed();
    if (index < size) return pointer{left, mirror, transpose};
    index -= size;
    return pointer{right, mirror, transpose};
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
  // assert(layer < nodes.size());
  // assert(pointer.index() < nodes[layer].size());
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

  auto descend = [&](auto layer, auto left, auto right) {
    const auto size = children(layer, left);
    const auto mirror = current.is_mirrored();
    const auto transpose = current.is_transposed();

    if (index < size) return pointer{left, mirror, transpose};
    index -= size;
    return pointer(right, mirror, transpose);
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

    for (const auto& node : nodes[layer]) {
      if (!node.left().empty())
        ++frequencies[node.left().index()];
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
      os << node << " (" << std::hash<detail::node>()(node) << "), ";
    }
    os << '\n';
  }

  os << "Leaves (" << nodes.front().size() << "): ";
  for (auto node : nodes.front())
    os << node.left().leaf() << ' ' << node.right().leaf() << " (" << std::hash<detail::node>()(node) << "), ";
  os << '\n';
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
    layer.emplace_back(index);
  } else {
    auto transform = created_node.transformations(canonical_node);
    layer.emplace_back(index, transform.first, transform.second);
  }
}