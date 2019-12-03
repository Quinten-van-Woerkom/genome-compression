#include "dna.h"

#include <array>
#include <bitset>
#include <cassert>
#include <cctype>
#include <iostream>
#include <string_view>

namespace fs = std::filesystem;

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
 *  Requires that <index> is smaller than <length>.
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

  nucleotide = std::toupper(nucleotide);

  switch (nucleotide) {
    case 'A': set_internal(index, false, false); break;
    case 'C': set_internal(index, true, false); break;
    case 'T': set_internal(index, false, true); break;
    case 'G': set_internal(index, true, true); break;
    // case '\n': case ' ': break; // We skip whitespace
    // default: {
    //   std::cerr << "Encountered unknown symbol, aborting...\n";
    //   exit(1);
    // }
    default: break; // We skip unknown symbols
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

auto read_genome(const fs::path path) -> std::vector<dna> {
  if (!fs::is_regular_file(path)) {
    std::cerr << "Non-existent path, aborting...\n";
    exit(1);
  }

  auto file = std::ifstream{path, std::ios::binary};
  auto buffer = std::array<char, dna::size()>{};
  auto result = std::vector<dna>{};

  if (!file.is_open()) {
    std::cerr << "Unable to open file, aborting...\n";
    exit(1);
  }

  while (true) {
    file.read(buffer.data(), buffer.size());
    if (file.eof()) break;
    result.emplace_back(std::string_view{buffer.data(), buffer.size()});
  }

  return result;
}