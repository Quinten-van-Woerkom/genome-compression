/**
 *  For now, the sequence is stored in terms of references a canonicalized,
 *  balanced binary tree.
 */

#pragma once

#include <iostream>
#include <random>
#include <string_view>
#include <tuple>

#include "robin_hood.h"

/******************************************************************************
 * FASTA nucleic acid codes
 *  Bit representations were chosen to allow certain bit twiddling hacks.
 *  Measurement shows that this does significantly improve performance.
 */
enum class nac : char {
  // Pairs of transpose base pairs
  A = 0b0001, T = 0b1000,
  C = 0b0010, G = 0b0100,
  R = 0b0011, Y = 0b1100,
  K = 0b0111, M = 0b1110,
  B = 0b0101, V = 0b1010,
  D = 0b1011, H = 0b1101,

  // These base pairs are their own transpose
  S = 0b0000, W = 0b1001,
  N = 0b0110, Indeterminate = 0b1111
};

/******************************************************************************
 * FASTA-compliant DNA strand
 *  Only Uracil is neglected, as it is not present in DNA; all other FASTA
 *  nucleic acid codes are supported.
 */
class dna {
public:
  dna() = default;
  dna(const std::string_view strand);
  dna(unsigned long long value) noexcept;

  static auto random(unsigned seed = 0) -> dna;
  static auto size() noexcept -> std::size_t { return length; }
  static auto size(std::size_t new_size) noexcept -> std::size_t { length = new_size; return length;}

  auto transposed() const noexcept -> dna;
  auto mirrored() const noexcept -> dna;
  auto inverted() const noexcept -> dna { return transposed().mirrored(); }
  auto invariant() const noexcept -> bool { return *this == mirrored(); }
  auto canonical() const noexcept -> std::tuple<dna, bool, bool, bool>;

  static auto bytes() noexcept -> std::size_t { return (size()+1)/2; }
  void serialize(std::ostream& os) const;
  static auto deserialize(std::istream& is) -> dna;

  auto code(std::size_t index) const -> nac;
  auto nucleotide(std::size_t index) const -> char;
  
  auto operator==(const dna& other) const noexcept -> bool { return nucleotides == other.nucleotides; }
  auto operator!=(const dna& other) const noexcept -> bool { return nucleotides != other.nucleotides; }
  auto operator<(const dna& other) const noexcept -> bool { return nucleotides < other.nucleotides; }

  operator std::uint64_t() const noexcept { return nucleotides; }
  auto to_ullong() const noexcept { return nucleotides; }

private:
  void set(std::size_t index, char nucleotide);
  void set(std::size_t index, nac code);

  std::uint64_t nucleotides;
  inline static std::size_t length = 12;   // Length of a single strand
};

auto operator<<(std::ostream& os, const dna& dna) -> std::ostream&;

namespace std {
  template<>
  struct hash<dna> {
    auto operator()(const dna& n) const noexcept -> std::size_t {
      return std::hash<std::uint64_t>()(n.to_ullong());
    }
  };
}