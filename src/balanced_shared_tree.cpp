/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 *  Specialised for balanced trees, to improve compression speed and ratio for
 *  those trees. This does mean that unbalanced trees are not supported.
 */

#include <limits>
#include <numeric>

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
 * Serializes the pointer to the given output stream in compressed format.
 * Starts with the 4 header bits and the most significant 4 bits, then
 * repeatedly stores the next most significant byte, until no more bytes
 * remain.
 */
void detail::pointer::serialize(std::ostream& os) const {
  int index = address_bits[segment]-4;
  std::uint8_t store = ((data >> index) & 0xf) | mirror << 4 | transpose << 5 | segment << 6;
  binary_write(os, store);
  for (index -= 8; index >= 0; index -= 8) {
    store = static_cast<std::uint8_t>(data >> index);
    binary_write(os, store);
  }
}

/**
 * Loads a pointer from an input stream.
 */
auto detail::pointer::deserialize(std::istream& is) -> pointer {
  auto result = pointer{};
  std::uint8_t loaded;
  binary_read(is, loaded);
  result.segment = (loaded >> 6) & 3u;
  result.transpose = (loaded >> 5) & 1u;
  result.mirror = (loaded >> 4) & 1u;
  
  auto index = address_bits[result.segment]-4;
  result.data = ((std::uint64_t)loaded & 0xf) << index;

  for (index -= 8; index >= 0; index -= 8) {
    binary_read(is, loaded);
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
auto balanced_shared_tree::access_leaf(pointer pointer) const -> dna {
  auto leaf = leaves[pointer.index()];
  if (pointer.is_mirrored()) leaf = leaf.mirrored();
  if (pointer.is_transposed()) leaf = leaf.transposed();
  return leaf;
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

/**
 * Indexing operator into the tree.
 * It is advised not to use this for iteration, as the induces overhead
 * compared to the implemented iterator is significant.
 */
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

  const auto& node = access_node(0, current);
  const auto mirror = current.is_mirrored();
  const auto transpose = current.is_transposed();

  const auto left = node.left();
  const auto right = node.right();
  if (index < !left.empty()) return access_leaf({left, mirror, transpose});
  else return access_leaf({right, mirror, transpose});
}

/**
 * Adds a leaf to the leaves layer.
 * Precondition: no similar leaves are already present in this layer.
 * Precondition: the leaf passed is of canonical variant.
 */
void balanced_shared_tree::emplace_leaf(dna leaf) {
  leaves.emplace_back(leaf);
}

/**
 * Adds a node to the specified layer.
 * Precondition: no similar nodes are already present in this layer. 
 */
void balanced_shared_tree::emplace_node(std::size_t layer, node node) {
  nodes[layer].emplace_back(node);
}

/**
 * Computes the number of times each pointer occurs in the nodes contained in
 * layer <layer>.
 * Precondition: <layer> is bigger than or equal to 0
 * Precondition: <layer> is smaller than the number of layers
 */
auto balanced_shared_tree::histogram(std::size_t layer) const -> std::vector<std::size_t> {
  assert(layer < nodes.size());
  const auto layer_size = (layer == 0) ? leaves.size() : nodes[layer-1].size();
  auto result = std::vector<std::size_t>(layer_size, 0);

  for (const auto& node : nodes[layer]) {
    if (auto left = node.left(); left) ++result[left.index()];
    if (auto right = node.right(); right) ++result[right.index()];
  }
  return result;
}

/**
 * Creates and stores a histogram for each layer in the tree, storing them as
 * lines in a .csv file.
 */
void balanced_shared_tree::store_histogram(std::filesystem::path path) const {
  auto file = std::ofstream{path};
  for (auto layer = 0u; layer < nodes.size(); ++layer) {
    auto frequencies = histogram(layer);
    std::sort(frequencies.begin(), frequencies.end(), std::greater<>());

    for (auto chunk : chunks(frequencies, 1000)) {
      for (const auto& frequency : chunk)
        file << frequency << ',';
      file << '\n';
    }
    file << '\n';
  }
}

/**
 * Sorts the pointers in each layer based on their relative reference count, to
 * reduce the pointers size required to refer to the most-referenced bits.
 * This further improves the effectiveness of pointer compression.
 */
void balanced_shared_tree::frequency_sort() {
  { // Leaf nodes are treated separately
    auto frequencies = histogram(0);
    auto old_indices = std::vector<std::size_t>(frequencies.size());
    std::iota(old_indices.begin(), old_indices.end(), 0);

    std::stable_sort(old_indices.begin(), old_indices.end(),
      [&](auto a, auto b) { return frequencies[a] > frequencies[b]; });

    auto new_leaves = std::vector<dna>(leaves.size());
    auto new_indices = std::vector<std::size_t>(old_indices.size());

    for (auto i = 0u; i < old_indices.size(); ++i)
      new_indices[old_indices[i]] = i;

    for (auto i = 0u; i < old_indices.size(); ++i)
      new_leaves[new_indices[i]] = leaves[i];
    leaves = std::move(new_leaves);

    auto update = [&](auto old) {
      if (old.empty()) return old;
      else return pointer{new_indices[old.index()], old.is_mirrored(), old.is_transposed()};
    };

    std::transform(nodes[0].begin(), nodes[0].end(), nodes[0].begin(),
      [&](auto old) { return node{update(old.left()), update(old.right())}; });
  }

  for (auto layer = 1u; layer < nodes.size(); ++layer) {
    auto frequencies = histogram(layer);
    auto old_indices = std::vector<std::size_t>(frequencies.size());
    std::iota(old_indices.begin(), old_indices.end(), 0);
    
    std::stable_sort(old_indices.begin(), old_indices.end(),
      [&](auto a, auto b) { return frequencies[a] > frequencies[b]; });
    
    auto new_nodes = nodes[layer-1];
    auto new_indices = old_indices;

    for (auto i = 0u; i < old_indices.size(); ++i)
      new_indices[old_indices[i]] = i;

    for (auto i = 0u; i < new_indices.size(); ++i)
      new_nodes[new_indices[i]] = nodes[layer-1][i];
    nodes[layer-1] = std::move(new_nodes);

    auto update = [&](auto old) {
      if (old.empty()) return old;
      else return pointer{new_indices[old.index()], old.is_mirrored(), old.is_transposed()};
    };

    std::transform(nodes[layer].begin(), nodes[layer].end(), nodes[layer].begin(),
      [&](auto old) { return node{update(old.left()), update(old.right())}; });
  }
}

/**
 * Computes the number of bytes required to store the compressed tree.
 */
auto balanced_shared_tree::bytes() const noexcept -> std::size_t {
  auto memory = root.bytes() + 8 + leaves.size()*sizeof(dna);

  for (const auto& layer : nodes) {
    memory += 8;  // Size of each layer is stored as 64 bits
    for (const auto& node : layer) memory += node.bytes();
  }
  return memory;
}

/**
 * Serializes the balanced tree to an output stream.
 * First stores the root, then all layers.
 * Each layer is stored as its length (as std::uint64_t), followed by all
 * separate nodes.
 */
void balanced_shared_tree::serialize(std::ostream& os) const {
  root.serialize(os);
  binary_write(os, leaves.size());
  for (const auto& leaf : leaves) leaf.serialize(os);

  for (const auto& layer : nodes) {
    binary_write(os, layer.size());
    for (const auto& node : layer) node.serialize(os);
  }
}

/**
 * Deserializes a balanced tree from an input stream.
 * Assumes it is stored starting with the root, followed by each layer, with
 * each layer stored as its size followed by the serialized nodes.
 */
auto balanced_shared_tree::deserialize(std::istream& is) -> balanced_shared_tree {
  auto result = balanced_shared_tree{};
  result.root = pointer::deserialize(is);
  std::uint64_t size;
  binary_read(is, size);
  for (auto i = 0u; i < size; ++i)
    result.leaves.emplace_back(dna::deserialize(is));

  while (true) {
    binary_read(is, size);
    if (!is) break;

    result.nodes.emplace_back();
    result.nodes.back().reserve(size);
    for (auto i = 0u; i < size; ++i)
      result.nodes.back().emplace_back(node::deserialize(is));
  }
  return result;
}

/**
 * Saves a balanced tree to a file in DAG format.
 */
void balanced_shared_tree::save(std::filesystem::path path) const {
  auto file = std::ofstream{path};
  serialize(file);
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
        auto mirror = top.is_mirrored();
        auto transpose = top.is_transposed();
        auto updated_pointer = pointer{next, mirror, transpose};
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
 * Checks if a leaf already exists in the tree, and if that is not the case,
 * inserts it into the map and into the tree dictionary.
 */
auto tree_constructor::emplace_leaf(dna leaf) -> pointer {
  auto [canonical, mirror, transpose] = leaf.canonical();
  auto insertion = leaves.emplace(canonical, parent.leaf_count());
  auto index = (*insertion.first).second;

  if (insertion.second) {
    parent.emplace_leaf(canonical);
  }
  return pointer{index, mirror, transpose};
}

/**
 * Emplaces leaves into the leaf map, if necessary, and adds a node referencing
 * them to the first non-leaf layer.
 */
auto tree_constructor::emplace_leaves(dna left, dna right) -> pointer {
  auto left_pointer = emplace_leaf(left);
  auto right_pointer = emplace_leaf(right);
  return emplace_node(0, left_pointer, right_pointer);
}

/**
 * Emplaces a single leaf in the map, and creates its parent node.
 * Used for leaves that have no neighbour on the right side.
 */
auto tree_constructor::emplace_leaves(dna last) -> pointer {
  auto pointer = emplace_leaf(last);
  return emplace_node(0, pointer);
}

/**
 * Constructs and emplaces a node inside the tree during its construction.
 * Returns a pointer to this node.
 */
auto tree_constructor::emplace_node(std::size_t layer, pointer left, pointer right) -> pointer
{
  auto created_node = node{left, right};
  auto insertion = nodes[layer].emplace(created_node, parent.node_count(layer));
  auto canonical_node = (*insertion.first).first;
  auto index = (*insertion.first).second;

  if (insertion.second) {
    parent.emplace_node(layer, created_node);
    return pointer{index, false, false};
  } else {
    auto transform = created_node.transformations(canonical_node);
    return pointer{index, transform.first, transform.second};
  }
}

/**
 * Reduces all gathered root nodes in order to fully reduce the tree.
 */
auto tree_constructor::reduce_roots() -> pointer {
  for (auto index = nodes.size(); roots.size() > 1; ++index)
    roots = reduce_nodes(roots, index);
  return roots.front();
}

/**
 * Reduces the current layer by constructing nodes out of the pointers it
 * contains, and emplacing those nodes in the correct layers and maps, if
 * necessary.
 */
auto tree_constructor::reduce_nodes(std::vector<pointer>& iterable, std::size_t index) -> std::vector<pointer> {
  auto layer = std::vector<pointer>{};
  layer.reserve(iterable.size()/2 + iterable.size()%2);

  if (parent.depth()-2 < index) {
    parent.add_layer();
    nodes.emplace_back();
  }

  foreach_pair(iterable,
    [&](auto left, auto right) { layer.emplace_back(emplace_node(index, left, right)); },
    [&](auto last) { layer.emplace_back(emplace_node(index, last)); }
  );

  return layer;
}