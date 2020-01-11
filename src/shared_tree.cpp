/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 */
 
 #include "shared_tree.h"

/******************************************************************************
 *  Annotated pointer to a node in the tree.
 */
auto pointer::data() const -> dna {
  assert(!empty() && "Trying to access null pointer.");
  assert(is_leaf() && "Trying to interpret a non-leaf node as a leaf.");
  return dna{information};
}

auto pointer::index() const -> std::size_t {
  assert(!empty() && "Trying to access null pointer.");
  assert(!is_leaf() && "Trying to interpret a leaf node as a non-leaf.");
  return information;
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
 *  Advances the iterator to the first node on the right of the current node.
 */
auto shared_tree::iterator::operator++() -> iterator& {
  while (!stack.empty()) {
    stack.pop_back();
    if (stack.back().is_leaf()) return *this;

    auto node = access(stack.back());
    stack.pop_back();
    if (auto right = node.right(); !right.empty()) stack.emplace_back(right);
    if (auto left = node.left(); !left.empty()) stack.emplace_back(left);
  }

  return *this;
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
  while (!current.is_leaf()) {
    const auto& node = access(current);
    const auto left_size = children(node.left());
    if (index < left_size) {
      current = node.left();
    } else {
      index -= left_size;
      current = node.right();
    }
  }
  return current.data();
}

/**
 *  Returns the node pointed to by the passed pointer.
 *  Accessing a nullptr results in undefined behaviour.
 */
auto shared_tree::access(pointer pointer) const -> const node& {
  const auto index = pointer.index();
  return nodes[index-1];
}