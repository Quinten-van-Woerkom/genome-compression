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

  node(pointer left, pointer right) : left{left}, right{right} {}
  node(pointer left) : left{left}, right{nullptr} {}

  // TODO: Add checks for empty children
  constexpr auto left_leaf() const -> Data { assert(left.is_leaf()); return left.get_leaf(); }
  constexpr auto right_leaf() const -> Data { assert(right.is_leaf()); return right.get_leaf(); }
  constexpr auto left_node() const -> const node& { assert(!left.is_leaf()); return left.get_node(); }
  constexpr auto right_node() const -> const node& { assert(!right.is_leaf()); return right.get_node(); }
  
  // Hash function
  constexpr auto hash() const noexcept -> std::size_t { return detail::hash(left, right); }

  constexpr auto size() const noexcept -> std::size_t { return left.size() + right.size(); }

  // For now, only mirror symmetry exists.
  // In the future, more similarities might be added.
  constexpr auto similarities(const node& other) const {
    assert(*this == other);
    if (left == other.right && right == other.left) return true;
    else return false;
  }

  // Accessor
  constexpr auto operator[](std::size_t index) const -> Data {
    if (index < left.size()) return left[index];
    else return right[index - left.size()];
  }

  constexpr auto operator==(const node& other) const noexcept {
    return left == other.left && right == other.right;
  }

  friend auto operator<<(std::ostream& os, const node& n) -> std::ostream& {
    return os << "node<" << n.left << ", " << n.right << ">";
  }

  /**
   *  Annotated pointer to another node in the tree.
   *  This node can be either a leaf or an internal node.
   *  Data is assumed to be approximately the same size as a pointer, so that
   *  it makes sense to store the data instead of a pointer.
   */
  class pointer {
  public:
    pointer() : pointer{nullptr} {}
    pointer(std::nullptr_t) : size_{0}, subnode{nullptr} {}  // Special construct to denote the absence of data
    pointer(const node& subnode) : size_{subnode.size()}, subnode{&subnode} {}
    pointer(Data data) : size_{1}, data{data} {}

    constexpr operator std::size_t() const noexcept { return detail::hash(as_unsigned, size_); }

    constexpr auto operator==(const pointer& other) const noexcept {
      return size_ == other.size_ && as_unsigned == other.as_unsigned;
    }

    // Accessor
    constexpr auto operator[](std::size_t index) const -> Data {
      assert(index < this->size());
      if (this->is_leaf()) return this->get_leaf();
      else return this->get_node()[index];
    }

    constexpr auto empty() const noexcept -> bool { return subnode == nullptr && !is_leaf(); }
    constexpr auto size() const noexcept -> std::size_t { return size_; }
    constexpr auto is_leaf() const noexcept -> bool { return size_ == 1; }
    // TODO: Add checks for emptiness
    constexpr auto get_leaf() const -> Data { assert(is_leaf() && "Trying to interpret a non-leaf node as a leaf."); return data; }
    constexpr auto get_node() const -> const node& { assert(!is_leaf() && "Trying to interpret a leaf node as a non-leaf."); return *subnode; }

    friend auto operator<<(std::ostream& os, const pointer& p) -> std::ostream& {
      if (p.empty()) return os << "empty";
      else if (p.is_leaf()) return os << "leaf: " << p.data;
      else return os << "node: " << p.subnode;
    }

  private:
    // TODO: Add annotations for similarity transforms
    std::size_t size_;
    
    union {
      const node* subnode;
      Data data;
      std::size_t as_unsigned;
    };
  };

private:
  pointer left, right;
};

namespace std {
  template<typename Data> struct hash<node<Data>> {
    auto operator()(const node<Data>& n) const noexcept -> std::size_t {
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
class shared_tree {
public:
  using node_type = node<Data>;
  using value_type = Data;
  using pointer = typename node_type::pointer;

  shared_tree(const std::vector<Data>& data);

  constexpr auto left_leaf() const -> Data { return root.get_node().left_leaf(); }
  constexpr auto right_leaf() const -> Data { return root.get_node().right_leaf(); }
  constexpr auto left_node() const -> const node_type& { return root.get_node().left_node(); }
  constexpr auto right_node() const -> const node_type& { return root.get_node().right_node(); }

  // Number of elements stored
  constexpr auto size() const -> std::size_t { return nodes.size(); }

  // Length of the data sequence
  constexpr auto length() const -> std::size_t { return root.size(); }

  constexpr auto operator[](std::size_t index) const -> Data {
    if (root.is_leaf()) return root.get_leaf();
    else return root.get_node()[index];
  }
  void print_unique() const {
    for (const auto& node : nodes)
        std::cout << node << '\n';
  }

private:
  pointer root;
  std::unordered_set<node_type> nodes;
};

/**
 *  Constructs a shared binary tree from a vector of data, using spatial
 *  subdivision for common subtree merging.
 */
template<typename Data>
shared_tree<Data>::shared_tree(const std::vector<Data>& data) {
    nodes.reserve(data.size()/2 + data.size()%2);
    using iter = typename std::unordered_set<node_type>::const_iterator;
    std::vector<iter> previous_layer;

    for (auto i = 0u; i < data.size() - 1; i += 2) {
      auto inserted = nodes.emplace(data[i], data[i+1]);
      previous_layer.emplace_back(inserted.first);
    }
    if (data.size() % 2) {
      auto inserted = nodes.emplace(data.back());
      previous_layer.emplace_back(inserted.first);
    }

    while (previous_layer.size() > 1) {
      std::vector<iter> next_layer;
      for (auto i = 0u; i < previous_layer.size() - 1; i += 2) {
        auto inserted = nodes.emplace(*previous_layer[i], *previous_layer[i+1]);
        next_layer.push_back(inserted.first);
      }
      if (previous_layer.size() % 2) {
        auto inserted = nodes.emplace(*previous_layer.back());
        next_layer.push_back(inserted.first);
      }
      previous_layer = std::move(next_layer);
    }

    root = pointer{*(previous_layer.front())};
  }