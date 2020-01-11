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

dna::dna(const std::string_view strand) {
  assert(strand.size() == length);
  for (auto i = 0u; i < length; ++i) {
    if (!valid_nac(strand[i])) std::cout << "Invalid strand: " << strand << '\n';
    set_nucleotide(i, strand[i]);
  }
}

dna::dna(unsigned long long value) noexcept : nucleotides{value} {}

/**
 *  Operator[] overload that redirects to return the nucleotide located at
 *  <index>, assuming that <index> is smaller than <length>.
 */
auto dna::operator[](std::size_t index) const -> char {
  return this->nucleotide(index);
}

/**
 *  Returns the nucleotide located at index <index>.
 *  Requires that <index> is smaller than <length>.
 */
auto dna::nucleotide(std::size_t index) const -> char {
  assert(index < length);
  auto nucleotide = static_cast<nac>(from_bits(nucleotides[4*index], nucleotides[4*index+1], nucleotides[4*index+2], nucleotides[4*index+3]));
  switch (nucleotide) {
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
 *  Internal helper function that sets the nucleotide located at <index> to
 *  <nucleotide>.
 */
void dna::set_nucleotide(std::size_t index, char nucleotide) {
  assert(index < length);

  auto set_internal = [&](auto index, auto nucleotide) {
    auto bits = to_bits(static_cast<unsigned char>(nucleotide));
    nucleotides.set(4*index, std::get<0>(bits));
    nucleotides.set(4*index+1, std::get<1>(bits));
    nucleotides.set(4*index+2, std::get<2>(bits));
    nucleotides.set(4*index+3, std::get<3>(bits));
  };

  nucleotide = std::toupper(nucleotide);

  switch (nucleotide) {
    case 'A': set_internal(index, nac::A); break;
    case 'C': set_internal(index, nac::C); break;
    case 'G': set_internal(index, nac::G); break;
    case 'T': set_internal(index, nac::T); break;
    case 'R': set_internal(index, nac::R); break;
    case 'Y': set_internal(index, nac::Y); break;
    case 'K': set_internal(index, nac::K); break;
    case 'M': set_internal(index, nac::M); break;
    case 'S': set_internal(index, nac::S); break;
    case 'W': set_internal(index, nac::W); break;
    case 'B': set_internal(index, nac::B); break;
    case 'D': set_internal(index, nac::D); break;
    case 'H': set_internal(index, nac::H); break;
    case 'V': set_internal(index, nac::V); break;
    case 'N': set_internal(index, nac::N); break;
    case '-': set_internal(index, nac::Indeterminate); break;
    default: {
      std::cerr << "Encountered unknown symbol: " << nucleotide << " (ASCII code " << static_cast<int>(nucleotide) << ")\n";
      exit(1);
    }
  }
}

/**
 *  Output stream operator that prints the nucleotides.
 */
auto operator<<(std::ostream& os, const dna& dna) -> std::ostream& {
  for (auto i = 0u; i < dna.size(); ++i) {
    os << dna.nucleotide(i);
  }
  return os;
}