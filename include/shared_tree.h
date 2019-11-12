/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 *  TODO: Implement symmetry
 */

#include <unordered_set>

namespace detail {
/// Hash function for any set of arguments.
/// Requires that the arguments be convertible to std::size_t
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
 *  Leaf node of a non-owning tree.
 *  Contains the stored data as well as a sentinel value to denote that it is a
 *  leaf node.
 */
template<typename Data>
class leaf {
public:

private:

}

/**
 *  Annotated pointer to another node in the tree.
 *  TODO: Allow similarity transforms, such as mirroring and inversion.
 */
template<typename Data>
class node {
public:

private:

}

/**
 *  Node in a non-owning tree. Children are either other nodes or leaf nodes.
 */
template<typename Data>
class internal_node {
public:

private:

}

/**
 *  Shared tree, where the manner of construction is not strictly given.
 *  Note that the tree itself is non-owning; rather, the nodes are stored
 *  elsewhere to allow for canonicalisation of nodes.
 */
template<typename Data>
class shared_tree {
public:


private:
    node* root;
}