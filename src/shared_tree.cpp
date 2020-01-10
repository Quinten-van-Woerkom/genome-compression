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
  return dna{data};
}

auto pointer::index() const -> std::size_t {
  assert(!empty() && "Trying to access null pointer.");
  assert(!is_leaf() && "Trying to interpret a leaf node as a non-leaf.");
  return data;
}


/******************************************************************************
 *  Genome, constructed as a balanced shared tree.
 */
auto shared_tree::children(pointer parent) const -> std::size_t {
  std::size_t total = 0llu;

  for (auto node : subtree(parent)) {
    ++total;
    (void)node;
  }

  return total;
}

/**
 *  Advances the iterator to the first node on the right of the current node.
 */
auto shared_tree::iterator::operator++() -> iterator& {
  while (!stack.empty()) {
    auto next = stack.back();
    stack.pop_back();

    if (next.is_leaf()) {
      current = next;
      // std::cout << "Found leaf: " << current << '\n';
      return *this;
    }

    // std::cout << "Accessing next: " << next << ' ' << access(next) << '\n';
    auto node = access(next);
    if (!node.right().empty()) {
      // std::cout << "Emplacing right: " << node.right() << '\n';
      stack.emplace_back(node.right());
    }
    if (!node.left().empty()) {
      // std::cout << "Emplacing left: " << node.left() << '\n';
      stack.emplace_back(node.left());
    }
  }

  current = nullptr;
  return *this;
}

/**
 *  Access function identical to the one in class shared_tree itself.
 */
auto shared_tree::iterator::access(pointer pointer) const -> const node& {
  const auto index = pointer.index();
  return (*nodes)[index-1];
}

auto shared_tree::operator[](std::size_t index) const -> dna {
  auto current = root;
  while (!current.is_leaf()) {
    const auto& node = nodes[current.index()-1];
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

auto shared_tree::access(pointer pointer) const -> const node& {
  const auto index = pointer.index();
  return nodes[index-1];
}