/**
 *  Buffered file reader implementation for single-FASTA DNA sequences.
 *  Allows one to read bigger DNA sequences without requiring them to be fully
 *  loaded in memory.
 */

#pragma once

#include <fstream>
#include <filesystem>
#include <iostream>
#include <string_view>
#include <thread>
#include <vector>

#include "dna.h"

class fasta_reader {
public:
  using value_type = dna;

  fasta_reader(std::filesystem::path path, std::size_t buffer_size = (1<<22));
  fasta_reader(const fasta_reader&) = delete;
  fasta_reader(fasta_reader&&) = default;

  ~fasta_reader() {
    if (background_loader.joinable()) background_loader.join();
  }

  auto eof() const -> bool { return end_of_file; }
  auto current_symbol() -> dna { return buffer[index]; }
  void next_symbol();
  void load_buffer();
  void swap_buffers();

  auto size() -> std::size_t;

  class iterator {
    public:
      iterator(fasta_reader& parent)
        : parent{parent} {}

      auto operator*() const -> dna { return parent.current_symbol(); }
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
  std::vector<dna> background_buffer;
  std::vector<dna> buffer;
  std::vector<char> char_buffer;
  std::ifstream file;
  std::filesystem::path path;
  std::size_t index;
  bool end_of_file = false;
  std::thread background_loader;
};

auto read_genome(const std::filesystem::path path) -> std::vector<dna>;