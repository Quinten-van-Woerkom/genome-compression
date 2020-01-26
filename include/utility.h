/**
 *  Implementation of convenience functions that are not strictly necessary
 *  but useful in cleaning up the code.
 */

#pragma once

#include <iostream>
#include <sstream>
#include <tuple>

/******************************************************************************
 *  Applies a functor to each consecutive pair. If the number of elements is
 *  odd, the last remaining entry is handled on its own.
 *  e.g. (1, 2, 3, 4, 5) -> (1, 2), (3, 4), (5)
 */
template<typename Iterable, typename BinaryFunc, typename UnaryFunc>
void foreach_pair(Iterable&& iterable, BinaryFunc binary_func, UnaryFunc unary_func) {
  auto begin = iterable.begin();
  auto end = iterable.end();

  while (begin != end) {
    auto&& left = *begin;
    ++begin;
    if (!(begin != end)) return unary_func(left);
    binary_func(left, *begin);
    ++begin;
  }
}

/******************************************************************************
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

/******************************************************************************
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


/******************************************************************************
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
      auto operator*() -> decltype(auto) { return *current; }
      auto operator*() const { return *current; }
      auto& operator++() { ++current; ++counter; return *this; }
      auto operator!=(const iterator& other) const { return counter < other.counter && current != other.current; }
      
      subiterator& current;
      std::size_t counter;
    };

    auto begin() { return iterator{begin_, 0}; }
    auto end() { return iterator{end_, chunk_size}; }
    auto size() const { return chunk_size; }

    subiterator& begin_, end_;
    std::size_t chunk_size;
  };

  struct iterator {
    auto operator*() { return chunk{begin, end, chunk_size}; }
    auto& operator++() { return *this; }
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


/******************************************************************************
 *  Constructs a range from an iterator pair.
 */
namespace detail {
template<typename Iterator1, typename Iterator2>
struct iterator_pair : std::pair<Iterator1, Iterator2> {
  using std::pair<Iterator1, Iterator2>::pair;

  auto begin() { return this->first; }
  auto end() { return this->second; }
  auto begin() const { return this->first; }
  auto end() const { return this->second; }
};
}

template<typename Iterator1, typename Iterator2>
auto iterator_pair(Iterator1 begin, Iterator2 end) {
  return detail::iterator_pair<Iterator1, Iterator2>{begin, end};
}


/******************************************************************************
 *  Parameter pack-supporting min.
 */
template<typename T>
constexpr auto variadic_min(T&& value) noexcept -> decltype(auto)
{
  return std::forward<T>(value);
}

template<typename T1, typename T2, typename... Ts>
constexpr auto variadic_min(T1&& value1, T2&& value2, Ts&&... values) noexcept -> decltype(auto)
{
  return value2 < value1
  ? variadic_min(value2, std::forward<Ts>(values)...)
  : variadic_min(value1, std::forward<Ts>(values)...);
}


/******************************************************************************
 *  File I/O helpers.
 *  Used for providing endianness-independent saving and loading.
 *  Data is stored in big-endian order.
 */
template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
void binary_write(std::ostream& os, T value, std::size_t bytes = sizeof(T)) {
  for (int i = bytes-1; i >= 0; --i) {
    unsigned char byte = (value >> (8*i)) & 0xff;
    os.put(byte);
  }
}

template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
void binary_read(std::istream& is, T& value, std::size_t bytes = sizeof(T)) {
  value = 0;
  for (int i = bytes-1; i >= 0; --i) {
    unsigned char byte;
    is.get(reinterpret_cast<char&>(byte));
    value |= ((T)byte << T(8*i));
  }
}


/******************************************************************************
 * Formats a size in bytes into the appropriate number of B, KB, MB, etc.
 */
template<typename T>
auto bytes_to_string(T bytes) {
  constexpr auto suffix = std::array{"B", "KB", "MB", "GB", "TB", "PB", "EB"};
  auto size = static_cast<double>(bytes);
  auto order = 0u;
  while (size >= 1000.0 && order < suffix.size()-1) {
    size /= 1000.0;
    ++order;
  }
  auto stream = std::stringstream{};
  stream.precision(3);
  stream << size << ' ' << suffix[order];
  return stream.str();
}


/******************************************************************************
 * Displays a progress bar filled to the given percentage.
 */
inline auto progress_bar(std::string_view name, unsigned current, unsigned end) {
  const auto percentage = double(current)/double(end);
  const auto bar_width = 60;
  const auto progress = unsigned(percentage*bar_width);
  auto output = std::stringstream{};
  output << '\r' << name << ": [";
  for (auto i = 0u; i < bar_width; ++i) {
    if (i < progress) output << '#';
    else output << ' ';
  }
  output << "] " << unsigned(percentage*100) << '%';
  return output.str();
}

inline auto spaces(unsigned length) {
  return std::string(length, ' ');
}