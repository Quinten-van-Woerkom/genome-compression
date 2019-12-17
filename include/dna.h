/**
 *  For now, the sequence is stored in terms of references a canonicalized,
 *  balanced binary tree.
 */

#pragma once

#include <array>
#include <bitset>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <string_view>
#include <vector>

namespace fs = std::filesystem;

// Supported nucleic acid codes
enum class nac : char {
  A, C, G, T, U, R, Y, K, M, S, W, B, D, H, V, N
};

// DNA strand of predetermined size
class dna {
  static constexpr std::size_t length = 1; // Length of a single strand
public:
  dna(const std::string_view strand);

  // Only for testing purposes
  static auto random(unsigned seed = 0) -> dna {
    dna result{};
    std::srand(seed);
    for (auto i = 0u; i < result.nucleotides.size(); ++i) {
      result.nucleotides.set(i, static_cast<bool>(std::rand()%2));
    }
    return result;
  }

  auto operator[](std::size_t index) const -> char;
  auto nucleotide(std::size_t index) const -> char;
  
  static constexpr auto size() noexcept -> std::size_t { return length; }

  friend auto operator<<(std::ostream& os, const dna& dna) -> std::ostream&;

  auto hash() const noexcept -> std::size_t { return std::hash<std::bitset<4*length>>{}(nucleotides); }
  auto operator==(const dna& other) const noexcept -> bool { return nucleotides == other.nucleotides; }
  auto operator!=(const dna& other) const noexcept -> bool { return nucleotides != other.nucleotides; }

private:
  dna() {}  // Does not initialize, useful only for testing purposes

  void set_nucleotide(std::size_t index, char nucleotide);

  std::bitset<4*length> nucleotides;
};

namespace std {
  template<>
  struct hash<dna> {
    auto operator()(const dna& n) const noexcept -> std::size_t {
      return n.hash();
    }
  };
}

auto read_genome(const fs::path path) -> std::vector<dna>;