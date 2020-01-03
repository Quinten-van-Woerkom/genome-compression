/**
 *  To be able to apply pointer compression, we require a self-contained heap
 *  so that we can actually sort the nodes based on frequency. This
 *  self-contained heap is implemented as a monotonic pool that supports
 *  ordinary indexing operations.
 */

#include <vector>

/**
 *  An arena stores a fixed-size buffer within which allocations can be done
 *  simply by incrementing a current pointer indicating the next available
 *  resource.
 */
template<typename T>
class arena {
public:
  using value_type = T;

  /**
   *  Instead of vector, we use manual new and delete, since the data in the
   *  buffer need not be constructed, that must be done by the user after
   *  allocation.
   */
  explicit arena(std::size_t size)
   : begin{new T[size]}, end{begin + size}, current{begin} {}

   ~arena() { delete[] begin; }

  /**
   *  Allocation is allowed to fail, as indicated by a nullptr return value,
   *  since a full arena is expected behaviour.
   */
  auto allocate(std::size_t size) noexcept -> T* {
    const auto remaining = end - current;
    if (size <= remaining) {
      const auto result = current;
      current += size;
      return result;
    } else return nullptr;
  }

  /**
   *  Since we expect a monotonic increase in memory resource, with
   *  simultaneous destruction, deallocation is a no-op. Memory is freed only
   *  when the arena itself is freed.
   */
  void deallocate(T* pointer, std::size_t size) const noexcept {}

  bool empty() const noexcept { return current == begin; }
  bool full() const noexcept { return current == end; }
  auto size() const noexcept { return end - begin; }
  bool contains(T* pointer) const noexcept { return pointer >= begin && pointer < end; }
  auto data() noexcept -> T* { return begin; }
  auto data() const noexcept -> const T* { return begin; }
  auto indexof(const T* pointer) const noexcept { return pointer - begin; }

  auto operator[](std::size_t index) const noexcept -> const T& { return begin[index]; }
  auto operator[](std::size_t index) noexcept -> T& { return begin[index]; }

private:
  T* begin;
  T* end;
  T* current;
};

/**
 *  Although we want some kind of illusion of contiguousness to allow for
 *  simple index-based accessors, we cannot guarantuee the existence of a
 *  single, contiguous region of memory that fits our data. As such, we
 *  implement a monotonic pool, where deallocation is a no-op.
 *  This slightly complicates the pointer arithmetic required to convert
 *  indices from and to pointers, but is more robust, and also likely to be
 *  faster when a lot of memory must be allocated.
 */
template<typename T>
class monotonic_pool {
public:
  using value_type = T;

  /**
   *  Constructs a monotonic pool, starting with the given amount and size of
   *  arenas. Note that the starting count and arena size must be at least one.
   */
  explicit monotonic_pool(std::size_t arena_size, std::size_t arena_count = 1)
    : arenas(arena_count, {arena_size}), arena_size{arena_size} {}

  /**
   *  Allocation is not allowed to fail for a monotonic pool, as opposed to an
   *  arena, as it is expected to grow indefinitely.
   */
  auto allocate(std::size_t size) -> T* {
    auto* result = arenas.back().allocate(size);
    if (result) return result;
    else {
      auto& a = arenas.emplace_back(arena_size);
      return a.allocate(size);
    }
  }

  /**
   *  Same as for the arenas, deallocation is a no-op.
   */
  void deallocate(T* pointer, std::size_t size) const noexcept {}

  bool empty() const { return arenas.front().empty(); }

  /**
   *  Computes the index of the object pointed to as if the data stored was
   *  stored contiguously. Throws a runtime error if the pointer is not
   *  contained in any of the arenas (not the best solution, but it works).
   */
  auto indexof(const T* pointer) const {
    for (auto i = 0u; i < arenas.size(); ++i)
      if (arenas[i].contains(pointer)) return arenas[i].indexof(pointer);
    throw std::runtime_error("Cannot find index of pointer not contained in memory pool");
  }

private:
  std::size_t arena_size;
  std::vector<arena<T>> arenas;
};