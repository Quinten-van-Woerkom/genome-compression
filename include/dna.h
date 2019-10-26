/**
 *  For now, the sequence is stored in terms of references to a context-free
 *  grammar, which is implemented as a canonicalized binary tree.
 */

#include <array>
#include <bitset>
#include <string_view>

// DNA strand of predetermined size
class dna {
private:
  static constexpr std::size_t length = 8; // Length of a single strand
public:
  dna(const std::string_view strand);

  auto operator[](std::size_t index) const -> char;
  auto nucleotide(std::size_t index) const -> char;
  
  constexpr auto size() const noexcept -> std::size_t { return length; }

  friend auto operator<<(std::ostream& os, const dna& dna) -> std::ostream&;

private:
  void set_nucleotide(std::size_t index, char nucleotide);

  std::bitset<2*length> nucleotides;
};