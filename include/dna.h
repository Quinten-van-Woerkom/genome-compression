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
#include <tuple>
#include <vector>

namespace fs = std::filesystem;

// Supported nucleic acid codes
enum class nac : char {
  A = 0, C, G, T, R, Y, K, M, S, W, B, D, H, V, N, Indeterminate
};

// DNA strand of predetermined size
class dna {
  static constexpr std::size_t length = 15; // Length of a single strand
public:
  dna(const std::string_view strand);
  dna(unsigned long long value) noexcept;

  static auto random(unsigned seed = 0) -> dna;
  static constexpr auto size() noexcept -> std::size_t { return length; }

  auto transposed() const noexcept -> dna;
  auto mirrored() const noexcept -> dna;
  auto inverted() const noexcept -> dna { return transposed().inverted(); }
  auto canonical() const noexcept -> std::tuple<dna, bool, bool>;
  auto to_ullong() const noexcept { return nucleotides.to_ullong(); }

  auto operator[](std::size_t index) const -> char;
  auto code(std::size_t index) const -> nac;
  auto nucleotide(std::size_t index) const -> char;
  
  auto operator==(const dna& other) const noexcept -> bool { return nucleotides == other.nucleotides; }
  auto operator!=(const dna& other) const noexcept -> bool { return nucleotides != other.nucleotides; }
  auto operator<(const dna& other) const noexcept -> bool { return nucleotides.to_ullong() < other.nucleotides.to_ullong(); }
  operator std::uint64_t() const noexcept { return nucleotides.to_ullong(); }

private:
  void set_nucleotide(std::size_t index, char nucleotide);
  void set_nucleotide(std::size_t index, nac code);

  std::bitset<4*length> nucleotides;
};

auto operator<<(std::ostream& os, const dna& dna) -> std::ostream&;

namespace std {
  template<>
  struct hash<dna> {
    auto operator()(const dna& n) const noexcept -> std::size_t {
      return std::hash<unsigned long long>{}(n.to_ullong());
    }
  };
}