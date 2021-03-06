/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 *  Specialised for balanced trees, to improve compression speed and ratio for
 *  those trees. This does mean that unbalanced trees are not supported.
 */

#include <cmath>
#include <future>
#include <limits>
#include <numeric>

#include "shared_tree.h"

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
  return 1ull << pointer::address_bits[segment];
}

/**
 * Helper function that returns the start of the address space associated with
 * a segment.
 */
constexpr auto address_start(std::size_t segment) noexcept {
  constexpr auto starts = std::array{
    0ull, address_space(0), address_space(0) + address_space(1),
    address_space(0) + address_space(1) + address_space(2)
  };
  return starts[segment];
}

/**
 * Helper function that determines the segment an index maps to.
 */
constexpr auto layer_segment(std::size_t index) noexcept {
  if (index < address_start(0b01)) return 0b00;
  else if (index < address_start(0b10)) return 0b01;
  else if (index < address_start(0b11)) return 0b10;
  else return 0b11;
}

/**
 * Helper function that determines the segment and offset an index belongs to.
 */
constexpr auto compress_pointer(std::size_t index) noexcept {
  if (index == 0x1fffffff) return std::pair{0b11, 0xfffffffull};
  const auto segment = layer_segment(index);
  return std::pair{segment, index - address_start(segment)};
}

/**
 * Helper function that determines the index that a segment-offset pair
 * corresponds with.
 */
constexpr auto decompress_pointer(std::size_t segment, std::size_t offset) noexcept {
  if (segment == 0b11 && offset == 0xfffffff) return 0x1fffffffull;
  return address_start(segment) + offset;
}

/**
 * Constructor from another pointer. Transformations can be applied to the
 * pointer, if necessary.
 * Note that a transposed nullptr must also be a nullptr, so that bits must
 * be set in this special case. This is the only case where a pointer is
 * invariant under transposition.
 */
pointer::pointer(const pointer& other, bool mirror_, bool transpose_) noexcept
: data{other.data},
  mirror{mirror_ != other.mirror && !other.invariant},
  transpose{transpose_ != other.transpose && other != nullptr},
  invariant{other.invariant} {}

/**
 * Constructor from an index with given similarity transforms.
 */
pointer::pointer(std::size_t index, bool mirror_, bool transpose_, bool invariant_) noexcept
: data{(std::uint32_t)index}, mirror{mirror_ && !invariant_}, transpose{transpose_}, invariant{invariant_} {}

/**
 * Constructs a null pointer, indicating an empty subtree.
 * Null pointers have all data bits set and are of maximum pointer size. Due to
 * their infrequent occurrence, this bigger size is not an issue, while it
 * significantly simplifies the indexing code.
 * Mirror and transpose are false for all null pointers and their
 * transformations.
 */
pointer::pointer(std::nullptr_t) noexcept
: data{0x1fffffff}, mirror{false}, transpose{false}, invariant{true} {}

/**
 * Returns the unsigned integer equivalent of the data stored in this
 * pointer. Note that the invariance bit is neglected.
 */
auto pointer::to_ulong() const noexcept -> unsigned long {
  return data
  | ((std::uint32_t)mirror << (address_bits.back() + 1))
  | ((std::uint32_t)transpose << (address_bits.back() + 2));
}


/**
 * Interprets the data as an index pointing to an inner node.
 */
auto pointer::index() const noexcept -> std::size_t {
  assert(!empty());
  return data;
}

/**
 * Returns the number of bytes required to store the pointer in memory, taking
 * into account the pointer compression applied on serialization.
 */
auto pointer::bytes() const noexcept -> std::size_t {
  const auto [segment, offset] = compress_pointer(data);
  return (4 + address_bits[segment])/8;
}

/**
 * Serializes the pointer to the given output stream in compressed format.
 * Starts with the 4 header bits and the most significant 4 bits, then
 * repeatedly stores the next most significant byte, until no more bytes
 * remain.
 */
void pointer::serialize(std::ostream& os) const {
  const auto [segment, offset] = compress_pointer(data);
  int index = address_bits[segment]-4;
  std::uint8_t store = offset >> index | mirror << 4 | transpose << 5 | segment << 6;
  binary_write(os, store);
  for (index -= 8; index >= 0; index -= 8) {
    store = static_cast<std::uint8_t>(offset >> index);
    binary_write(os, store);
  }
}

/**
 * Loads a pointer from an input stream.
 */
auto pointer::deserialize(std::istream& is) -> pointer {
  std::uint8_t loaded;
  binary_read(is, loaded);
  auto segment = (loaded >> 6) & 0x3;
  bool transpose = (loaded >> 5) & 0x1;
  bool mirror = (loaded >> 4) & 0x1;
  auto index = address_bits[segment]-4;
  auto offset = ((std::uint64_t)loaded & 0xf) << index;

  for (index -= 8; index >= 0; index -= 8) {
    binary_read(is, loaded);
    offset |= ((std::uint64_t)loaded << index);
  }

  auto data = decompress_pointer(segment, offset);
  return pointer{data, mirror, transpose, false};
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
bool node::operator==(const node& other) const noexcept {
  return children == other.children;
}

/**
 * Serializes the node to the given output stream.
 * First, the left pointer is stored, and then the right pointer.
 */
void node::serialize(std::ostream& os) const {
  children[0].serialize(os);
  children[1].serialize(os);
}

/**
 * Deserializes a node from an input stream.
 * The node is assumed to be stored in compressed format.
 */
auto node::deserialize(std::istream& is) -> node {
  auto left = pointer::deserialize(is);
  auto right = pointer::deserialize(is);
  return node{left, right};
}


/******************************************************************************
 * class shared_tree:
 *  Shared binary tree class that exploits structural properties of balanced
 *  trees to store its directed acyclic graph representation more effectively.
 */
/**
 * Constructs a shared_tree from a FASTA formatted file.
 */
shared_tree::shared_tree(fasta_reader file, bool verbose) {
  auto constructor = tree_constructor{*this};
  root = constructor.reduce(file, verbose);
}

shared_tree::shared_tree(std::vector<dna>& data, bool verbose) {
  auto constructor = tree_constructor{*this};
  root = constructor.reduce(data, verbose);
}

/**
 * Returns the total number of nodes in the tree.
 * Leafs are excluded.
 */
auto shared_tree::node_count() const -> std::size_t {
  auto sum = 0ull;
  for (const auto& layer : nodes) sum += layer.size();
  return sum;
}

/**
 * Accesses the first or second leaf located in the node pointed to by
 * <pointer>, depending on <index>.
 */
auto shared_tree::access_leaf(pointer pointer) const -> dna {
  auto leaf = leaves[pointer.index()];
  if (pointer.is_mirrored()) leaf = leaf.mirrored();
  if (pointer.is_transposed()) leaf = leaf.transposed();
  return leaf;
}

/**
 * Accesses the node in layer <layer> pointed to by <pointer>.
 * Returned by value since since the nodes must be immutable anyway.
 */
auto shared_tree::access_node(std::size_t layer, pointer pointer) const -> node {
  return nodes[layer][pointer.index()];
}

/**
 * Returns the number of children contained in the subtree referenced by
 * <pointer> in <layer>. To determine this, we just traverse the tree without
 * considering mirroring or transposition, as those do not alter the number of
 * children a node has.
 */
auto shared_tree::children(std::size_t layer, pointer pointer) const -> std::size_t {
  if (pointer.empty()) return 0;
  const auto node = access_node(layer, pointer);
  const auto left = node.left();
  const auto right = node.right();
  if (layer == 0) return !left.empty() + !right.empty();
  else return children(layer-1, left) + children(layer-1, right); 
}

/**
 * Indexing operator into the tree.
 * It is advised not to use this for iteration, as the induced overhead
 * compared to the implemented iterator is significant.
 * Precondition: index < width
 * Precondition: no nullptrs within the tree, only at the right edge
 */
auto shared_tree::operator[](std::uint64_t index) const -> dna {
  auto current = root;

  auto descend = [&](auto layer, auto left, auto right) {
    const auto size = (1u << (layer));
    const auto mirror = current.is_mirrored();
    const auto transpose = current.is_transposed();

    if (index < size) return pointer{left, mirror, transpose};
    index -= size;
    return pointer(right, mirror, transpose);
  };

  for (int layer = nodes.size()-1; layer >= 0; --layer) {
    const auto& node = access_node(layer, current);
    const auto left = node.left();
    const auto right = node.right();

    if (current.is_mirrored()) current = descend(layer, right, left);
    else current = descend(layer, left, right);
  }

  return access_leaf(current);
}

/**
 * Adds a leaf to the leaves layer.
 * Precondition: no similar leaves are already present in this layer.
 * Precondition: the leaf passed is of canonical variant.
 */
void shared_tree::emplace_leaf(dna leaf) {
  leaves.emplace_back(leaf);
}

/**
 * Adds a node to the specified layer.
 * Precondition: no similar nodes are already present in this layer. 
 */
void shared_tree::emplace_node(std::size_t layer, node node) {
  nodes[layer].emplace_back(node);
}

/**
 * Computes the number of times each pointer occurs in the nodes contained in
 * layer <layer>.
 * Precondition: <layer> is bigger than or equal to 0
 * Precondition: <layer> is smaller than the number of layers
 */
auto shared_tree::histogram(std::size_t layer) const -> std::vector<std::size_t> {
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
void shared_tree::store_histogram(std::filesystem::path path) const {
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
 * Helper functions in layer sorting.
 */
auto iota(std::size_t size) {
  auto vector = std::vector<std::size_t>(size);
  std::iota(vector.begin(), vector.end(), 0);
  return vector;
}

/**
 * Inverts a vector of indices, so that the value i located at j
 * becomes the value of j located at i.
 */
auto invert_indices(const std::vector<std::size_t>& indices) {
  auto inverted = indices;
  for (auto i = 0u; i < indices.size(); ++i)
    inverted[indices[i]] = i;
  return inverted;
};

/**
 * Returns the child layer, but then reordered so that child i is now located at
 * index indices[i].
 */
template<typename T>
auto reorder_layer(const std::vector<T>& children, const std::vector<std::size_t>& indices) {
  auto reordered = children;
  for (auto i = 0u; i < indices.size(); ++i)
    reordered[indices[i]] = children[i];
  return reordered;
};

/**
 * Rewires all nodes to point to the correct children according to the child
 * reshuffling as indicated by indices.
 */
void shared_tree::rewire_nodes(std::size_t layer, const std::vector<std::size_t>& indices) {
  // Rewires a pointer to point to the same child, but then sorted.
  auto rewire_pointer = [](auto old, const auto& indices) {
    if (old.empty()) return old;
    const auto index = indices[old.index()];
    const auto mirror = old.is_mirrored();
    const auto transpose = old.is_transposed();
    const auto invariant = old.is_invariant();
    return pointer{index, mirror, transpose, invariant};
  };

  // Rewires a node to point to its children, but sorted.
  auto rewire_node = [&](auto old, const auto& indices) {
    auto left = rewire_pointer(old.left(), indices);
    auto right = rewire_pointer(old.right(), indices);
    return node{left, right};
  };

  std::transform(nodes[layer].begin(), nodes[layer].end(), nodes[layer].begin(),
    [&](auto old) { return rewire_node(old, indices); });
}

/**
 * Sorts the leaves based on frequency of reference by the parent layer.
 * Also rewires the parent nodes to match this shuffle.
 */
void shared_tree::sort_leaves() {
  auto frequencies = histogram(0);
  auto indices = iota(frequencies.size());

  std::stable_sort(indices.begin(), indices.end(),
    [&](auto a, auto b) { return frequencies[a] > frequencies[b]; });


  indices = invert_indices(indices);
  leaves = reorder_layer(leaves, indices);
  rewire_nodes(0, indices);
}

/**
 * Sorts the layer based on frequency of reference by its parent layer.
 * Also rewires those parent nodes to match this shuffle.
 */
void shared_tree::sort_nodes(std::size_t layer) {
  auto frequencies = histogram(layer+1);
  auto indices = iota(frequencies.size());

  std::stable_sort(indices.begin(), indices.end(),
    [&](auto a, auto b) { return frequencies[a] > frequencies[b]; });

  indices = invert_indices(indices);
  nodes[layer] = reorder_layer(nodes[layer], indices);
  rewire_nodes(layer+1, indices);
}

/**
 * Sorts the pointers in each layer based on their relative reference count, to
 * reduce the pointers size required to refer to the most-referenced bits.
 * This further improves the effectiveness of pointer compression.
 */
void shared_tree::sort_tree(bool verbose) {
  if (verbose)
    std::cout << progress_bar("Sorting nodes", 0, 1) << std::flush;
  std::vector<std::future<void>> futures;

  // Checks if a future is ready
  auto ready = [](const auto& future) {
    if (!future.valid()) return false;
    return future.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
  };

  futures.emplace_back(std::async(&shared_tree::sort_leaves, this));
  for (auto layer = 1u; layer < nodes.size()-1; layer += 2)
    futures.emplace_back(std::async(&shared_tree::sort_nodes, this, layer));

  auto count = 0;
  for (auto& future : futures) {
    if (ready(future)) {
      future.get();
      ++count;
      if (verbose)
        std::cout << progress_bar("Sorting nodes", count, 2*futures.size()) << std::flush;
    }
  }

  futures.clear();
  for (auto layer = 0u; layer < nodes.size()-1; layer += 2)
    futures.emplace_back(std::async(&shared_tree::sort_nodes, this, layer));

  for (auto& future : futures) {
    if (ready(future)) {
      future.get();
      ++count;
      if (verbose)
        std::cout << progress_bar("Sorting nodes", count, 2*futures.size()) << std::flush;
    }
  }

  if (verbose)
    std::cout << "\rSorting nodes: done." << spaces(100) << '\n';
}

/**
 * Computes the number of bytes required to store the compressed tree.
 */
auto shared_tree::bytes() const noexcept -> std::size_t {
  auto memory = root.bytes() + 8 + leaves.size()*dna::bytes();

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
void shared_tree::serialize(std::ostream& os) const {
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
auto shared_tree::deserialize(std::istream& is) -> shared_tree {
  auto result = shared_tree{};
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
void shared_tree::save(std::filesystem::path path) const {
  auto file = std::ofstream{path};
  serialize(file);
}


/******************************************************************************
 * class shared_tree::iterator:
 *  Iterator over the tree.
 */
shared_tree::iterator::iterator(shared_tree& parent, std::size_t layer, pointer root) : parent{parent} {
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
auto shared_tree::iterator::operator*() const noexcept -> dna {
  auto top = stack.back().current;
  return parent.access_leaf(top);
}

/**
 *  Advances the iterator to the first leaf on the right of the current leaf.
 *  To achieve this, pops the last-found leaf off the stack and recurses on the
 *  remaining stack members until a leaf is found.
 */
auto shared_tree::iterator::operator++() -> iterator& {
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
void shared_tree::iterator::next_leaf() {
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
tree_constructor::tree_constructor(shared_tree& parent)
: parent{parent} {
  nodes.reserve(64);
}

/**
 * Checks if a leaf already exists in the tree, and if that is not the case,
 * inserts it into the map and into the tree dictionary.
 */
auto tree_constructor::emplace_leaf(dna leaf) -> pointer {
  const auto [canonical, mirror, transpose, invariant] = leaf.canonical();
  const auto insertion = leaves.emplace(canonical, parent.leaf_count());
  const auto index = (*insertion.first).second;

  if (insertion.second) parent.emplace_leaf(canonical);
  return pointer{index, mirror, transpose, invariant};
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
  const auto created_node = node{left, right};
  const auto [canonical_node, mirror, transpose] = created_node.canonical();
  const auto insertion = nodes[layer].emplace(canonical_node, parent.node_count(layer));
  if (insertion.second) parent.emplace_node(layer, canonical_node);
  const auto index = (*insertion.first).second;

  const auto invariant = (left == right.mirrored());
  return pointer{index, mirror, transpose, invariant};
}

/**
 * Reduces all gathered root nodes in order to fully reduce the tree.
 */
auto tree_constructor::reduce_roots(bool verbose) -> pointer {
  const auto size = log2(roots.size());
  auto i = 0;
  for (auto index = nodes.size(); roots.size() > 1; ++index, ++i) {
    roots = reduce_nodes(roots, index);
    if (verbose)
      std::cout << progress_bar("Combining subtrees", i, size) << std::flush;
  }

  if (verbose)
    std::cout << "\rCombining subtrees: done." << spaces(100) << '\n';

  return roots.front();
}

/**
 * Reduces the current layer by constructing nodes out of the pointers it
 * contains, and emplacing those nodes in the correct layers and maps, if
 * necessary.
 */
auto tree_constructor::reduce_nodes(const std::vector<pointer>& iterable, std::size_t index) -> std::vector<pointer> {
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

/**
 * Reduces data read from a file into segments, each of which is fully reduced.
 * The resulting subtree roots are then accumulated into a single top layer
 * which is also reduced to complete the tree.
 */
auto tree_constructor::reduce(fasta_reader& file, bool verbose) -> pointer {
  std::vector<dna> buffer;
  auto current_buffer = 0;
  const auto approximate_buffer_count = file.buffers();
  while (file.read_into(buffer)) {
    reduce_segment(buffer);

    if (verbose) {
      ++current_buffer;
      std::cout << progress_bar("Constructing subtrees", current_buffer, approximate_buffer_count) << std::flush;
    }
  }

  if (verbose)
    std::cout << "\rConstructing subtrees: done." << spaces(100) << '\n';

  return reduce_roots(verbose);
}

/**
 * Reduces the data by dividing it into segments, each of which is fully
 * reduced. The resulting tree roots are then also reduced to obtain the
 * final tree representation.
 */
auto tree_constructor::reduce(const std::vector<dna>& data, bool verbose) -> pointer {
  constexpr auto subtree_depth = 25;
  constexpr auto subtree_width = (1u<<subtree_depth);

  auto current_subtree = 0;
  const auto subtrees = data.size() / subtree_width;

  for (auto segment : chunks(data, subtree_width)) {
    reduce_segment(segment);
    
    if (verbose) {
      ++current_subtree;
      std::cout << progress_bar("Constructing subtrees", current_subtree, subtrees) << std::flush;
    }
  }

  if (verbose)
    std::cout << "\rConstructing subtrees: done." << spaces(100) << '\n';
  
  return reduce_roots(verbose);
}