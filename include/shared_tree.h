/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 */

#include <unordered_set>

namespace detail {
/**
 *  Hash function for an arbitrary set of arguments.
 *  Requires only that each argument be convertible to std::size_t.
 */
template <typename... Args> auto hash(Args &&... args) noexcept -> std::size_t {
  if constexpr (sizeof...(args) == 0) return 0;
  else {
    return hash_impl(args...);
  }
}

template <typename Arg, typename... Args>
constexpr auto hash_impl(Arg&& arg, Args&& ... args) noexcept -> std::size_t {
  constexpr auto scalar = (1 << (sizeof...(args) + 1)) + 1;
  return scalar * arg + hash_impl(args...);
}

constexpr auto hash_impl() noexcept -> std::size_t { 
    return 0;
}
}

/**
 *  Annotated pointer to another node in the tree.
 *  This node can be either a leaf or an internal node.
 *  Data is assumed to be approximately the same size as a pointer.
 */
template<typename Data>
class node_pointer {
public:
  // TODO: Add accessor functions


  operator std::size_t() const noexcept { return as_unsigned; }

private:
  // TODO: Add annotations for similarity transforms
  bool is_leaf : 1;
  union {
    node* subnode;
    Data data;
    std::size_t as_unsigned;
  }
}

/**
 *  Node in a non-owning tree. Children are either other nodes or leaf nodes,
 *  containing the actual data stored.
 */
template<typename Data>
class node {
public:

private:
  node_pointer<Data> left, right;
}

/**
 *  Shared tree, where the manner of construction is not strictly given.
 *  Note that the tree itself is non-owning; rather, the nodes must be stored
 *  elsewhere to allow for canonicalisation of nodes.
 */
template<typename Data>
class shared_tree {
public:


private:
  node_pointer root_node;
}