/**
 *  Shared tree, i.e. a tree compressed to a directed acyclic graph through
 *  common subtree merging.
 */
 
 #include "shared_tree.h"

/******************************************************************************
 *  Annotated pointer to a node in the tree.
 */
auto pointer::size() const noexcept -> std::size_t {
  if (empty()) return 0;
  if (is_leaf()) return 1;
  return subnode->size();
}

auto pointer::get_leaf() const -> const dna& {
  assert(is_leaf() && "Trying to interpret a non-leaf node as a leaf.");
  return data;
}

auto pointer::get_node() const -> const node& {
  assert(!is_leaf() && "Trying to interpret a leaf node as a non-leaf.");
  return *subnode;
}


/******************************************************************************
 *  Genome, constructed as a balanced shared tree.
 */
auto shared_tree::operator[](std::size_t index) const -> const dna& {
  const auto* current = &root;
  while (!current->is_leaf()) {
    const auto& node = current->get_node();
    if (index < node.left().size()) {
      current = &node.left();
    } else {
      index -= node.left().size();
      current = &node.right();
    }
  }
  return current->get_leaf();
}