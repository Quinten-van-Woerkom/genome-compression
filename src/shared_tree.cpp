/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 */
 
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
 *  Constructs an iterator over the shared tree.
 *  A nullptr argument for the root indicates an end iterator.
 */
shared_tree::iterator::iterator(const std::vector<node>& nodes, pointer root) : nodes{nodes} {
  if (root != nullptr) {
    stack.emplace_back(root);
    next_leaf();
  }
}

/**
 *  Advances the iterator to the first node on the right of the current node.
 *  To achieve this, pops the last-found leaf off the stack and recurses on the
 *  remaining stack members.
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
    auto top = stack.back();
    if (top.is_leaf()) return;

    auto node = access(top);
    stack.pop_back();
    if (auto right = node.right(); !right.empty()) stack.emplace_back(right);
    if (auto left = node.left(); !left.empty()) stack.emplace_back(left);
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
  return current.leaf();
}

/**
 *  Returns the node pointed to by the passed pointer.
 *  Accessing a nullptr results in undefined behaviour.
 */
auto shared_tree::access(pointer pointer) const -> const node& {
  const auto index = pointer.index();
  return nodes[index-1];
}