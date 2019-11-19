/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 */

#include <cassert>
#include <unordered_set>
#include <vector>

namespace detail {
/**
 *  Hash function for an arbitrary set of arguments.
 *  Requires only that each argument be convertible to std::size_t.
 */
constexpr auto hash_impl() noexcept -> std::size_t {
  return 0;
}

template <typename Arg, typename... Args>
constexpr auto hash_impl(const Arg& arg, const Args& ... args) noexcept -> std::size_t {
  constexpr auto scalar = (1 << (sizeof...(args) + 1)) + 1;
  return scalar * arg + hash_impl(args...);
}

template <typename... Args> auto hash(const Args&... args) noexcept -> std::size_t {
  if constexpr (sizeof...(args) == 0) return 0;
  else {
    return hash_impl(args...);
  }
}
}

/**
 *  Node in a non-owning tree. Children are either other nodes or leaf nodes,
 *  containing the actual data stored.
 */
template<typename Data>
class node {
public:
  class pointer;

  node(node& left, node& right) : left{left}, right{right} {}
  node(Data left, Data right) : left{left}, right{right} {}
  node(node& left) : left{left}, right{nullptr} {}
  node(Data left) : left{left}, right{nullptr} {}

  // TODO: Add checks for empty children
  constexpr auto left_leaf() -> Data { assert(left.is_leaf()); return left.get_leaf(); }
  constexpr auto right_leaf() -> Data { assert(right.is_leaf()); return right.get_leaf(); }
  constexpr auto left_node() -> node& { assert(!left.is_leaf()); return left.get_node(); }
  constexpr auto right_node() -> node& { assert(!right.is_leaf()); return right.get_node(); }
  
  // Hash function
  constexpr auto hash() const noexcept -> std::size_t { return detail::hash(left, right); }

  constexpr auto operator==(const node& other) const noexcept {
    return left == other.left && right == other.right;
  }

  /**
   *  Annotated pointer to another node in the tree.
   *  This node can be either a leaf or an internal node.
   *  Data is assumed to be approximately the same size as a pointer, so that
   *  it makes sense to store the data instead of a pointer.
   */
  class pointer {
  public:
    pointer(std::nullptr_t) : leaf_node{false}, subnode{nullptr} {}  // Special construct to denote the absence of data
    pointer(node& subnode) : leaf_node{false}, subnode{&subnode} {}
    pointer(Data data) : leaf_node{true}, data{data} {}

    constexpr operator std::size_t() const noexcept { return as_unsigned; }

    constexpr auto operator==(const pointer& other) const noexcept {
      return leaf_node == other.leaf_node && as_unsigned == other.as_unsigned;
    }

    constexpr auto empty() const noexcept -> bool { return subnode == nullptr && !is_leaf(); }
    constexpr auto is_leaf() const noexcept -> bool { return leaf_node; }
    // TODO: Add checks for emptiness
    constexpr auto get_leaf() -> Data { assert(is_leaf() && "Trying to interpret a non-leaf node as a leaf."); return data; }
    constexpr auto get_node() -> node& { assert(!is_leaf() && "Trying to interpret a leaf node as a non-leaf."); return *subnode; }

  private:
    // TODO: Add annotations for similarity transforms
    bool leaf_node : 1;
    union {
      node* subnode;
      Data data;
      std::size_t as_unsigned;
    };
  };

private:
  pointer left, right;
};

namespace std {
  template<typename Data> struct hash<node<Data>> {
    std::size_t operator()(const node<Data>& n) const noexcept {
      return n.hash();
    }
  };
}

/**
 *  Shared tree, where the manner of construction is not strictly given.
 *  Note that the tree itself is non-owning; rather, the nodes must be stored
 *  elsewhere to allow for canonicalisation of nodes.
 *  Precondition: each node is assumed to exist only once.
 */
template<typename Data>
class shared_binary_tree {
public:
  using node_type = node<Data>;
  using pointer = typename node_type::pointer;

  shared_binary_tree(pointer root, std::unordered_set<node_type> nodes)
    : root{root}, nodes{nodes} {}

  static auto construct_from(const std::vector<Data>& data) -> shared_binary_tree {
    std::unordered_set<node_type> nodes;
    using iter = typename std::unordered_set<node_type>::iterator;
    std::vector<iter> previous_layer;

    for (auto i = 0u; i < data.size() - 1; i += 2) {
      auto new_node = nodes.emplace(data[i], data[i+1]).first;
      previous_layer.push_back(new_node);
    }
    if (data.size() % 2 != 0) {
      auto new_node = nodes.emplace(data.back()).first;
      previous_layer.push_back(new_node);
    }

    // while (previous_layer.size() > 1) {
    //   std::vector<iter> next_layer;
    //   for (auto i = 0u; i < data.size() - 1; i += 2) {
    //     auto new_node = nodes.emplace(*previous_layer[i], *previous_layer[i+1]).first;
    //     next_layer.push_back(new_node);
    //   }
    //   if (data.size() % 2 != 0) {
    //     auto new_node = nodes.emplace(*previous_layer.back()).first;
    //     next_layer.push_back(new_node);
    //   }
    //   previous_layer = next_layer;
    // }

    auto root = pointer{*(previous_layer.front())};
    return shared_binary_tree{root, nodes};
  }

private:
  pointer root;
  std::unordered_set<node_type> nodes;
};