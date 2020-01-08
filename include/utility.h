/**
 *  Implementation of convenience functions that are not strictly necessary
 *  but useful in cleaning up the code.
 */

#pragma once

#include <iostream>
#include <tuple>

/**
 *  Range that iterates over two ranges simultaneously.
 *  TODO: Seems to be bugged for non-copyable objects, to be fixed.
 */
namespace detail {
template<typename Iterable1, typename Iterable2>
class zip_class {
public:
  zip_class(Iterable1 it1, Iterable2 it2)
    : it1{std::forward<Iterable1>(it1)}, it2{std::forward<Iterable2>(it2)} {}

  template<typename Iterator1, typename Iterator2>
  class iterator {
  public:
    iterator(Iterator1 it1, Iterator2 it2) : it1{it1}, it2{it2} {}

    auto operator*() { return std::forward_as_tuple(*it1, *it2); }
    auto& operator++() { ++it1; ++it2; return *this; }
    auto operator++(int) { auto temp = *this; ++*this; return temp; }
    auto& operator--() { --it1; --it2; return *this; }
    auto operator--(int) { auto temp = *this; --*this; return temp; }
    auto operator!=(const iterator& other) { return it1 != other.it1 && it2 != other.it2; } // Stop iteration at first end

  private:
    Iterator1 it1;
    Iterator2 it2;
  };

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
auto zip(Iterable1&& it1, Iterable2&& it2) -> detail::zip_class<Iterable1&&, Iterable2&&> {
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
  struct iterator {
    iterator(Iterator left, Iterator right) : left{left}, right{right} {}

    auto operator*() { return std::forward_as_tuple(*left, *right); }
    auto& operator++() { ++left; ++left; ++right; ++right; return *this; }
    auto operator++(int) { auto temp = *this; ++left; ++left; ++right; ++right; return temp; }
    auto& operator--() { --left; --left; --right; --right; return *this; }
    auto operator--(int) { auto temp = *this; --left; --left; --right; --right; return temp; }

    // We stop iteration when our iterator matches the end's left iterator.
    auto operator!=(const iterator& other) {
      return left != other.left && right != other.left;
    }

    // Whether or not the iterators overlap exactly, i.e. it is an exact match.
    // This would indicate that the range contained an even number of elements
    // and could therefore be visited in full.
    bool exactly_equals(const iterator& other) {
      return !(left != other.left || right != other.right);
    }

    Iterator left, right;
  };

  auto begin() { return iterator{iterable.begin(), ++iterable.begin()}; }
  auto end() { return iterator{iterable.end(), ++iterable.end()}; }
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

template<typename Iterable, typename BinaryFunc, typename UnaryFunc>
void pairwise_apply(Iterable&& iterable, BinaryFunc binary_func, UnaryFunc unary_func) {
  auto&& range = pairwise(std::forward<Iterable>(iterable));
  auto begin = range.begin();
  auto end = range.end();
  for (; begin != end; ++begin)
    binary_func(*begin.left, *begin.right);
  
  if (!begin.exactly_equals(end))
    unary_func(*begin.left);
}

/**
 *  Hash function for an arbitrary set of arguments.
 *  Requires only that each argument be convertible to std::size_t.
 */
namespace detail {
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
 *  Compose unsigned integers from bits, and decompose them in their bits.
 *  Bits are order from least significant to most significant bit.
 */
constexpr auto from_bits(bool bit) noexcept -> unsigned long long {
  return bit;
}

template<typename... T>
constexpr auto from_bits(bool low_bit, T... bits) noexcept -> unsigned long long {
  return low_bit | (from_bits(bits...) << 1);
}

template<typename T>
constexpr auto to_bits(T value) noexcept -> std::array<bool, 8*sizeof(T)> {
  std::array<bool, 8*sizeof(T)> result = {};
  for (auto i = 0u; i < 8*sizeof(T); ++i) {
    result[i] = value & 1;
    value = value >> 1;
  }
  return result;
}


/**
 *  Range adaptor that returns the range divided into chunks.
 */
namespace detail {
template<typename Iterable>
class chunks_class {
public:
  chunks_class(Iterable iterable, std::size_t chunk_size) : iterable{iterable}, chunk_size{chunk_size} {}

  using subiterator = decltype(std::begin(std::declval<Iterable>()));
  using value_type = decltype(*std::declval<subiterator>());

  struct chunk {
    struct iterator {
      auto operator*() { return std::forward<value_type>(*current); }
      auto operator*() const { return *current; }
      auto& operator++() { ++current; ++counter; return *this; }
      auto operator!=(const iterator& other) const { return counter < other.counter && current != other.current; }
      
      subiterator current;
      std::size_t counter;
    };

    auto begin() const { return iterator{current, 0}; }
    auto end() const { return iterator{end_, chunk_size}; }
    auto size() const { return chunk_size; }

    subiterator current, end_;
    std::size_t chunk_size;
  };

  struct iterator {
    auto operator*() const { return chunk{begin, end, chunk_size}; }
    auto& operator++() { for (auto i = 0u; i < chunk_size && begin != end; ++i) ++begin; return *this; }
    auto operator!=(const iterator&) const { return begin != end; }

    subiterator begin, end;
    std::size_t chunk_size;
  };

  auto begin() { return iterator{std::begin(iterable), std::end(iterable), chunk_size}; }
  auto end() { return iterator{std::begin(iterable), std::end(iterable), chunk_size}; }


private:
  Iterable&& iterable;
  std::size_t chunk_size;
};
}

template<typename Iterable>
auto chunks(Iterable&& iterable, std::size_t chunk_size) -> detail::chunks_class<Iterable> {
  return {std::forward<Iterable>(iterable), chunk_size};
}