/**
 *  Buffered file reader implementation for single-FASTA DNA sequences.
 *  Allows one to read bigger DNA sequences without requiring them to be fully
 *  loaded in memory.
 */

#pragma once

#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <vector>

class fasta_reader {
public:
  fasta_reader(std::filesystem::path file, std::size_t buffer_size = 32768, std::size_t step_size = 1);
  fasta_reader(const fasta_reader&) = delete;
  fasta_reader(fasta_reader&&) = default;

  auto eof() const -> bool { return file.eof() && (index >= current_buffer.size()); }
  void next_symbol();
  void next_character();
  void skip_comment();
  auto current_symbol() -> std::string_view;
  auto current_character() const -> char { return current_buffer[index]; }

  // Merely an upper bound, as headers are also considered.
  // As each character is a single base pair, the number of base pairs is
  // equivalent to the size of the file.
  auto size() -> std::size_t {
    auto pos = file.tellg();
    file.seekg(0);
    file.ignore(std::numeric_limits<std::streamsize>::max());
    auto length = file.gcount();
    file.clear();
    file.seekg(pos);
    return length;
  }

  class iterator {
    public:
      iterator(fasta_reader& parent)
        : parent{parent} {}

      auto operator*() const -> std::string_view { return parent.current_symbol(); }
      auto operator==(const iterator&) const -> bool { return parent.eof(); } // Hackish implementation but it works
      auto operator!=(const iterator&) const -> bool { return !parent.eof(); }

      auto operator++() {
        parent.next_symbol();
        return *this;
      }

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
  bool verbose;
};