#include "dna.h"

#include <array>
#include <bitset>
#include <cassert>
#include <cctype>
#include <iostream>
#include <string_view>
#include <tuple>

#include "fasta_reader.h"
#include "utility.h"

namespace fs = std::filesystem;

bool valid_nac(char code) {
  code = std::toupper(code);
  return code == 'A' || code == 'C' || code == 'G' ||code == 'T'
    || code == 'R' || code == 'Y' || code == 'K' || code == 'M'
    || code == 'S' || code == 'W' || code == 'B' || code == 'D'
    || code == 'H' || code == 'V' || code == 'N' || code == '-';
}

auto to_nac(char nucleotide) -> nac {
  auto code = std::toupper(nucleotide);
  switch (code) {
    case 'A': return nac::A;
    case 'C': return nac::C;
    case 'G': return nac::G;
    case 'T': return nac::T;
    case 'R': return nac::R;
    case 'Y': return nac::Y;
    case 'K': return nac::K;
    case 'M': return nac::M;
    case 'S': return nac::S;
    case 'W': return nac::W;
    case 'B': return nac::B;
    case 'D': return nac::D;
    case 'H': return nac::H;
    case 'V': return nac::V;
    case 'N': return nac::N;
    case '-': return nac::Indeterminate;
    default : {
      std::cerr << "Encountered unknown symbol: " << code << " (ASCII code " << static_cast<int>(code) << ")\n";
      exit(1);
    }
  }
}

constexpr auto from_nac(nac code) noexcept -> char {
  switch (code) {
    case nac::A: return 'A';
    case nac::C: return 'C';
    case nac::G: return 'G';
    case nac::T: return 'T';
    case nac::R: return 'R';
    case nac::Y: return 'Y';
    case nac::K: return 'K';
    case nac::M: return 'M';
    case nac::S: return 'S';
    case nac::W: return 'W';
    case nac::B: return 'B';
    case nac::D: return 'D';
    case nac::H: return 'H';
    case nac::V: return 'V';
    case nac::N: return 'N';
    case nac::Indeterminate: return '-';
    default: {
      std::cerr << "Trying to decipher invalid FASTA symbol\n";
      exit(1);
    }
  }
}

constexpr auto transpose(nac code) noexcept -> nac {
  switch (code) {
    case nac::A: return nac::T;
    case nac::T: return nac::A;
    case nac::C: return nac::G;
    case nac::G: return nac::C;
    case nac::R: return nac::Y;
    case nac::Y: return nac::R;
    case nac::K: return nac::M;
    case nac::M: return nac::K;
    case nac::S: return nac::S;
    case nac::W: return nac::W;
    case nac::B: return nac::V;
    case nac::V: return nac::B;
    case nac::D: return nac::H;
    case nac::H: return nac::D;
    case nac::N: return nac::N;
    case nac::Indeterminate: return nac::Indeterminate;
    default: {
      std::cerr << "Trying to transpose an unknown nucleic acid code\n";
      exit(1);
    }
  }
}

/**
 * Constructors.
 */
dna::dna(const std::string_view strand) {
  assert(strand.size() == length);
  for (auto i = 0u; i < length; ++i)
    set_nucleotide(i, strand[i]);
}

dna::dna(unsigned long long value) noexcept : nucleotides{value} {}

/**
 *  Returns a random-initialised DNA strand. Used for testing purposes.
 */
auto dna::random(unsigned seed) -> dna {
  std::srand(seed);
  return dna{static_cast<unsigned long long>(rand() | ((std::uint64_t)rand() << 32))};
}

/**
 * Returns a transposed version of the DNA strand.
 */
auto dna::transposed() const noexcept -> dna {
  auto result = *this;
  for (auto i = 0u; i < dna::size(); ++i)
    result.set_nucleotide(i, transpose(code(i)));
  return result;
}

/**
 * Returns a mirrored version of the DNA strand.
 */
auto dna::mirrored() const noexcept -> dna {
  auto result = *this;
  for (auto i = 0u; i < dna::size(); ++i)
    result.set_nucleotide(i, code(dna::size() - i - 1));
  return result;
}

/**
 * Returns the canonical node representation of this DNA sequence, as well as
 * the transformations necessary to obtain it from the current representation.
 * The canonical version is determined to be the one with the lowest bit
 * representation.
 * The booleans returned indicate the requirement of mirroring and/or
 * transformation to transform from the current to the canonical
 * representation.
 * The third boolean represents whether or not a strand is invariant under
 * mirroring.
 */
auto dna::canonical() const noexcept -> std::tuple<dna, bool, bool, bool> {
  const auto is_invariant = invariant();
  auto current = std::tuple{*this, false, false, is_invariant};
  const auto transpose = std::tuple{transposed(), false, true, is_invariant};
  const auto mirror = std::tuple{mirrored(), true, false, is_invariant};
  const auto both = std::tuple{transposed().mirrored(), true, true, is_invariant};

  if (transpose < current) current = transpose;
  if (mirror < current) current = mirror;
  if (both < current) current = both;
  return current;
}

/**
 * Serializes the DNA strand into an output stream.
 * Little-endian storage format is used.
 */
void dna::serialize(std::ostream& os) const {
  binary_write(os, nucleotides.to_ullong());
}

/**
 * Deserializes the DNA strand from an input stream.
 * Little-endian storage format is used.
 */
auto dna::deserialize(std::istream& is) -> dna {
  std::uint64_t value;
  binary_read(is, value);
  return dna{value};
}

/**
 *  Operator[] overload that redirects to return the nucleotide located at
 *  <index>, assuming that <index> is smaller than <length>.
 */
auto dna::operator[](std::size_t index) const -> char {
  return this->nucleotide(index);
}

/**
 * Returns the nucleic acid code of the nucleotide located at <index>
 * Requires that <index> is smaller than <length>.
 */
auto dna::code(std::size_t index) const -> nac {
  assert(index < length);
  return static_cast<nac>(from_bits(nucleotides[4*index], nucleotides[4*index+1], nucleotides[4*index+2], nucleotides[4*index+3]));
}

/**
 *  Returns the nucleotide located at index <index>.
 *  Requires that <index> is smaller than <length>.
 */
auto dna::nucleotide(std::size_t index) const -> char {
  assert(index < length);
  auto nac = code(index);
  return from_nac(nac);
}

/**
 *  Internal helper function that sets the nucleotide located at <index> to
 *  <nucleotide>.
 */
void dna::set_nucleotide(std::size_t index, char nucleotide) {
  auto code = to_nac(nucleotide);
  set_nucleotide(index, code);
}

void dna::set_nucleotide(std::size_t index, nac code) {
  assert(index < length);
  auto bits = to_bits(static_cast<unsigned char>(code));
  nucleotides.set(4*index, std::get<0>(bits));
  nucleotides.set(4*index+1, std::get<1>(bits));
  nucleotides.set(4*index+2, std::get<2>(bits));
  nucleotides.set(4*index+3, std::get<3>(bits));
}

/**
 *  Output stream operator that prints the nucleotides.
 */
auto operator<<(std::ostream& os, const dna& dna) -> std::ostream& {
  for (auto i = 0u; i < dna.size(); ++i)
    os << dna.nucleotide(i);
  os << " (" << dna.invariant() << ')';
  return os;
}