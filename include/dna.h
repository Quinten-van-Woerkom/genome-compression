/**
 *  For now, the sequence is stored in terms of references to a context-free
 *  grammar, which is implemented as a canonicalized binary tree.
 */

#include <array>
#include <bitset>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string_view>
#include <vector>

namespace fs = std::filesystem;

// DNA strand of predetermined size
class dna {
  static constexpr std::size_t length = 8; // Length of a single strand
public:
  dna(const std::string_view strand);

  auto operator[](std::size_t index) const -> char;
  auto nucleotide(std::size_t index) const -> char;
  
  static constexpr auto size() noexcept -> std::size_t { return length; }

  friend auto operator<<(std::ostream& os, const dna& dna) -> std::ostream&;

  auto hash() const noexcept -> std::size_t { return nucleotides.to_ulong(); }
  auto operator==(const dna& other) const noexcept -> bool { return nucleotides == other.nucleotides; }

private:
  void set_nucleotide(std::size_t index, char nucleotide);

  std::bitset<2*length> nucleotides;
};

namespace std {
  template<>
  struct hash<dna> {
    std::size_t operator()(const dna& n) const noexcept {
      return n.hash();
    }
  };
}

auto read_genome(const fs::path path) -> std::vector<dna>;