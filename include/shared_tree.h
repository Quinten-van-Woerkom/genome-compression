/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 */

#include <cassert>
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
 *  Node in a non-owning tree. Children are either other nodes or leaf nodes,
 *  containing the actual data stored.
 */
template<typename Data>
class node {
private:
  class pointer;

public:
  node(const node& left, const node& right) : left{left}, right{right} {}
  node(Data left, Data right) : left{left}, right{right} {}

  constexpr auto left_leaf() const -> Data { assert(left.is_leaf()); return left.get_leaf(); }
  constexpr auto right_leaf() const -> Data { assert(right.is_leaf()); return right.get_leaf(); }
  constexpr auto left_node() const -> node& { assert(!left.is_leaf()); return left.get_node(); }
  constexpr auto right_node() const -> node& { assert(!right.is_leaf()); return right.get_node(); }
  
private:
  /**
   *  Annotated pointer to another node in the tree.
   *  This node can be either a leaf or an internal node.
   *  Data is assumed to be approximately the same size as a pointer.
   */
  class pointer {
  public:
    pointer(const node& subnode) : leaf_node{false}, subnode{subnode} {}
    pointer(Data data) : leaf_node{true}, data{data} {}

    // TODO: Add accessor functions
    constexpr operator std::size_t() const noexcept { return as_unsigned; }

    constexpr auto is_leaf() const noexcept -> bool { return leaf_node; }
    constexpr auto get_leaf() const -> Data { assert(is_leaf() && "Trying to interpret a non-leaf node as a leaf."); return data; }
    constexpr auto get_node() const -> node& { assert(!is_leaf() && "Trying to interpret a leaf node as a non-leaf."); return *subnode; }

  private:
    // TODO: Add annotations for similarity transforms
    bool leaf_node : 1;
    union {
      node* subnode;
      Data data;
      std::size_t as_unsigned;
    };
  };

  pointer left, right;
};

/**
 *  Shared tree, where the manner of construction is not strictly given.
 *  Note that the tree itself is non-owning; rather, the nodes must be stored
 *  elsewhere to allow for canonicalisation of nodes.
 *  Precondition: each node is assumed to exist only once.
 */
template<typename Data>
class shared_tree {
public:


private:
  using pointer = typename node<Data>::pointer;
  pointer root_node;
};