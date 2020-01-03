/**
 *  Buffered file reader implementation for single-FASTA DNA sequences.
 *  Allows one to read bigger DNA sequences without requiring them to be fully
 *  loaded in memory.
 */

#pragma once

#include <cassert>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string_view>
#include <vector>

class fasta_reader {
public:
  fasta_reader(std::filesystem::path path, std::size_t step_size = 1, std::size_t buffer_size = 32768);
  fasta_reader(const fasta_reader&) = delete;
  fasta_reader(fasta_reader&&) = default;

  auto eof() const -> bool { return file.eof() && index >= (buffer.size()-1); }
  auto current_symbol() -> std::string_view { return {buffer.data() + index, step_size}; }
  void next_symbol();
  void load_buffer();

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

      auto& operator++() {
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
  std::vector<char> buffer;
  std::ifstream file;
  std::size_t index;
  std::size_t step_size;
  bool wrapped_comment = false;
};