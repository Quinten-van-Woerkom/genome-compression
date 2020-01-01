/**
 *  Buffered file reader implementation for single-FASTA DNA sequences.
 *  Allows one to read bigger DNA sequences without requiring them to be fully
 *  loaded in memory.
 */

#include "file_reader.h"

#include <iostream>

namespace fs = std::filesystem;

fasta_reader::fasta_reader(fs::path file, std::size_t buffer_size, std::size_t step_size)
  : subbuffer(step_size), current_buffer(buffer_size), previous_buffer(buffer_size),
    file{file, std::ios::binary}, swapped{false}, index{0} {
  if (!this->file.is_open()) {
    std::cerr << "Unable to open file, aborting...\n";
    exit(1);
  }

  this->file.read(current_buffer.data(), current_buffer.size());
  this->next_symbol();
}

void fasta_reader::next_symbol() {
  swapped = false;
  index += subbuffer.size();
  for (auto i = 0u; i < subbuffer.size(); ++i) {
    next_character();
    if (current_character() == '>') skip_comment();
  }
}

void fasta_reader::next_character() {
  if (index >= current_buffer.size()) {
    swapped = true;
    std::swap(previous_buffer, current_buffer);
    file.read(current_buffer.data(), current_buffer.size());
    if (file.eof()) current_buffer.resize(file.gcount());
    index -= previous_buffer.size();
  }
}

void fasta_reader::skip_comment() {
  do {
    next_character();
  } while (current_character() != '\n' || current_character() == '>');
}

auto fasta_reader::current_symbol() -> std::string_view {
  if (!swapped) return std::string_view{&current_buffer[index-subbuffer.size()], subbuffer.size()};
  std::copy(&previous_buffer[previous_buffer.size() - subbuffer.size() + index], previous_buffer.data() + previous_buffer.size(), subbuffer.data());
  std::copy(&current_buffer[0], &current_buffer[index], &subbuffer[subbuffer.size() - index]);
  return std::string_view{subbuffer.data(), subbuffer.size()};
}