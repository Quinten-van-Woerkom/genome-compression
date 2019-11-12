#include "dna.h"

#include <array>
#include <bitset>
#include <cassert>
#include <iostream>
#include <string_view>

dna::dna(const std::string_view strand) {
  assert(strand.size() == length);
  for (auto i = 0u; i < length; ++i) {
    set_nucleotide(i, strand[i]);
  }
}

/**
 *  Operator[] overload that redirects to return the nucleotide located at
 *  <index>, assuming that <index> is smaller than <length>.
 */
auto dna::operator[](std::size_t index) const -> char {
  return this->nucleotide(index);
}

/**
 *  Returns the nucleotide located at index <index>.
 *  Assumes that <index> is smaller than <length>.
 */
auto dna::nucleotide(std::size_t index) const -> char {
  assert(index < length);
  auto nucleotide = nucleotides.test(2*index) | (nucleotides.test(2*index+1) << 1);
  switch (nucleotide) {
    case 0b00: return 'A';
    case 0b01: return 'C';
    case 0b10: return 'T';
    case 0b11: return 'G';
    default: return 'N';
  }
}

/**
 *  Internal helper function that sets the nucleotide located at <index> to
 *  <nucleotide>.
 */
void dna::set_nucleotide(std::size_t index, char nucleotide) {
  assert(index < length);

  auto set_internal = [&](auto index, bool low, bool high) {
    nucleotides.set(2*index, low);
    nucleotides.set(2*index+1, high);
    assert(nucleotides.test(2*index) == low);
    assert(nucleotides.test(2*index+1) == high);
  };

  switch (nucleotide) {
    case 'A': set_internal(index, false, false); break;
    case 'C': set_internal(index, true, false); break;
    case 'T': set_internal(index, false, true); break;
    case 'G': set_internal(index, true, true); break;
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