#include "dna.h"

#include <cassert>
#include <cctype>
#include <iostream>
#include <string_view>
#include <tuple>

#include "fasta_reader.h"
#include "utility.h"

namespace fs = std::filesystem;

/******************************************************************************
 * Helper functions for FASTA parsing.
 */
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

constexpr auto from_nac(nac code) -> char {
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

/**
 * Constructors.
 */
dna::dna(const std::string_view strand) : nucleotides{0} {
  assert(strand.size() == length);
  for (auto i = 0u; i < length; ++i) {
    set(i, strand[i]);
  }
}

dna::dna(unsigned long long value) noexcept : nucleotides{value} {}

/**
 *  Returns a random-initialised DNA strand. Used for testing purposes.
 */
auto dna::random(unsigned seed) -> dna {
  std::srand(seed);
  auto random = static_cast<unsigned long long>(rand() | ((std::uint64_t)rand() << 32));
  auto mask = (1u << dna::size()) - 1;
  return dna{random & mask};
}

/**
 * Returns a transposed version of the DNA strand.
 * The DNA representation is exploited to allow simple nibble inversion to be
 * equivalent to transposition. Rather unnecessary but measurements show it
 * does affect performance.
 */
auto dna::transposed() const noexcept -> dna {
  auto v = nucleotides;
  // swap odd and even bits
  v = ((v >> 1) & 0x5555555555555555) | ((v & 0x5555555555555555) << 1);
  // swap consecutive pairs
  v = ((v >> 2) & 0x3333333333333333) | ((v & 0x3333333333333333) << 2);
  return dna{v};
}

/**
 * Returns a mirrored version of the DNA strand.
 */
auto dna::mirrored() const noexcept -> dna {
  dna result{};
  for (auto i = 0u; i < length; ++i)
    result.set(i, code(length - i - 1));
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
 * mirroring. To determine this, we use the fact that all similar nodes are
 * invariant if any one is.
 */
auto dna::canonical() const noexcept -> std::tuple<dna, bool, bool, bool> {
  const auto is_invariant = invariant();
  const auto current = std::tuple{*this, false, false, is_invariant};
  const auto transpose = std::tuple{transposed(), false, true, is_invariant};
  const auto mirror = std::tuple{mirrored(), true, false, is_invariant};
  const auto invert = std::tuple{inverted(), true, true, is_invariant};

  return variadic_min(current, transpose, mirror, invert);
}

/**
 * Serializes the DNA strand into an output stream.
 * Big-endian storage format is used.
 */
void dna::serialize(std::ostream& os) const {
  binary_write(os, nucleotides, (dna::size()+1)/2);
}

/**
 * Deserializes the DNA strand from an input stream.
 * Big-endian storage format is used.
 */
auto dna::deserialize(std::istream& is) -> dna {
  std::uint64_t value = 0;
  binary_read(is, value, (dna::size() + 1)/2);
  return dna{value};
}

/**
 * Returns the nucleic acid code of the nucleotide located at <index>
 * Requires that <index> is smaller than <length>.
 */
auto dna::code(std::size_t index) const -> nac {
  assert(index < length);
  const auto offset = 4*index;
  return static_cast<nac>((nucleotides >> offset) & 0xf);
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
void dna::set(std::size_t index, char nucleotide) {
  auto code = to_nac(nucleotide);
  set(index, code);
}

void dna::set(std::size_t index, nac code) {
  assert(index < length);
  const auto offset = 4*index;
  nucleotides &= ~(0xfull << offset);
  nucleotides |= static_cast<std::uint64_t>(code) << offset;
}

/**
 *  Output stream operator that prints the nucleotides.
 */
auto operator<<(std::ostream& os, const dna& strand) -> std::ostream& {
  for (auto i = 0u; i < strand.size(); ++i)
    os << strand.nucleotide(i);
  return os;
}