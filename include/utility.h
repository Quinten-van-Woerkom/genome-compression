/**
 *  Implementation of convenience functions that are not strictly necessary
 *  but useful in cleaning up the code.
 */

#pragma once

#include <iostream>

/**
 *  Range that iterates over two ranges simultaneously.
 */
namespace detail {
template<typename Iterable1, typename Iterable2>
class zip_class {
public:
  zip_class(Iterable1 it1, Iterable2 it2)
    : it1{it1}, it2{it2} {}

  template<typename Iterator1, typename Iterator2>
  class iterator {
  public:
    iterator(Iterator1 it1, Iterator2 it2) : it1{it1}, it2{it2} {}

    auto operator*() { return std::pair{*it1, *it2}; }
    auto& operator++() { ++it1; ++it2; return *this; }
    auto operator++(int) { auto temp = *this; ++it1; ++it2; return temp; }
    auto& operator--() { --it1; --it2; return *this; }
    auto operator--(int) { auto temp = *this; --it1; --it2; return temp; }
    auto operator!=(const iterator& other) { return it1 != other.it1 && it2 != other.it2; } // Stop iteration at first end

  private:
    Iterator1 it1;
    Iterator2 it2;
  };

  // template<typename Iterator1, typename Iterator2>
  // class const_iterator {
  // public:
  //   const_iterator(Iterator1 it1, Iterator2 it2) : it1{it1}, it2{it2} {}

  //   auto operator*() { return std::pair{*it1, *it2}; }
  //   auto& operator++() { ++it1; ++it2; return *this; }
  //   auto operator++(int) { auto temp = *this; ++it1; ++it2; return temp; }
  //   auto& operator--() { --it1; --it2; return *this; }
  //   auto operator--(int) { auto temp = *this; --it1; --it2; return temp; }
  //   auto operator!=(const const_iterator& other) { return it1 != other.it1 && it2 != other.it2; } // Stop iteration at first end

  // private:
  //   Iterator1 it1;
  //   Iterator2 it2;
  // };

  auto begin() { return iterator{it1.begin(), it2.begin()}; }
  auto end() { return iterator{it1.end(), it2.end()}; }
  auto cbegin() const { return iterator{it1.cbegin(), it2.cbegin()}; }
  auto cend() const { return iterator{it1.cend(), it2.cend()}; }

private:
  Iterable1 it1;
  Iterable2 it2;
};
}

template<typename Iterable1, typename Iterable2>
auto zip(Iterable1&& it1, Iterable2&& it2) -> detail::zip_class<Iterable1, Iterable2> {
  return {std::forward<Iterable1>(it1), std::forward<Iterable2>(it2)};
}


/**
 *  Range that iterates over a sequence in pairs.
 */
namespace detail {
template<typename Iterable>
class pairwise_class {
public:
  pairwise_class(Iterable iterable) : iterable{iterable} {}

  template<typename Iterator>
  class iterator {
  public:
    iterator(Iterator left, Iterator right) : left{left}, right{right} {}

    auto operator*() { return std::pair{*left, *right}; }
    auto& operator++() { ++left; ++left; ++right; ++right; return *this; }
    auto operator++(int) { auto temp = *this; ++left; ++left; ++right; ++right; return temp; }
    auto& operator--() { --left; --left; --right; --right; return *this; }
    auto operator--(int) { auto temp = *this; --left; --left; --right; --right; return temp; }

    // We stop iteration when our iterator matches the end's left iterator.
    auto operator!=(const iterator& other) { return left != other.left && right != other.left; }

  private:
    Iterator left, right;
  };

  // template<typename Iterator>
  // class const_iterator {
  // public:
  //   const_iterator(Iterator left, Iterator right) : left{left}, right{right} {}

  //   auto operator*() { return std::pair{*left, *right}; }
  //   auto& operator++() { ++left; ++left; ++right; ++right; return *this; }
  //   auto operator++(int) { auto temp = *this; ++left; ++left; ++right; ++right; return temp; }
  //   auto& operator--() { --left; --left; --right; --right; return *this; }
  //   auto operator--(int) { auto temp = *this; --left; --left; --right; --right; return temp; }
  //   auto operator!=(const iterator& other) { return left != other.left && left != other.right; }

  // private:
  //   Iterator left, right;
  // };

  auto begin() { return iterator{iterable.begin(), std::next(iterable.begin())}; }
  auto end() { return iterator{iterable.end(), std::next(iterable.end())}; }
  auto cbegin() { return iterator{iterable.cbegin(), ++iterable.cbegin()}; }
  auto cend() { return iterator{iterable.cend(), ++iterable.cend()}; }

private:
  Iterable iterable;
};
}

template<typename Iterable>
auto pairwise(Iterable&& iterable) -> detail::pairwise_class<Iterable> {
  return {std::forward<Iterable>(iterable)};
}