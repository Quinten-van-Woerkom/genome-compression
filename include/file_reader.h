/**
 *  Buffered file reader implementation for single-FASTA DNA sequences.
 *  Allows one to read bigger DNA sequences without requiring them to be fully
 *  loaded in memory.
 */

#include <fstream>
#include <filesystem>
#include <string_view>
#include <vector>

class fasta_reader {
public:
  fasta_reader(std::filesystem::path file, std::size_t buffer_size = 32768, std::size_t step_size = 1);

  auto eof() const -> bool { return file.eof() && (index >= current_buffer.size()); }
  void next();
  void step();
  auto current() -> std::string_view;
  auto character() const -> char { return current_buffer[index]; }

  class iterator {
    public:
      constexpr iterator(fasta_reader& parent) : parent{parent} {}
      auto operator*() const -> std::string_view { return parent.current(); }
      auto operator==(const iterator&) const -> bool { return parent.eof(); } // Hackish implementation but it works
      auto operator!=(const iterator&) const -> bool { return !parent.eof(); }
      auto operator++() { parent.next(); return *this; }
    private:
      fasta_reader& parent;
  };

  using const_iterator = iterator;

  auto begin() { return iterator{*this}; }
  auto end() { return iterator{*this}; }

private:
  std::vector<char> subbuffer;
  std::vector<char> current_buffer;
  std::vector<char> previous_buffer;
  std::ifstream file;
  bool swapped;
  std::size_t index;
};