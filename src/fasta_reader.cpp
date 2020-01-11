/**
 *  Buffered file reader implementation for single-FASTA DNA sequences.
 *  Allows one to read bigger DNA sequences without requiring them to be fully
 *  loaded in memory.
 */

#include "fasta_reader.h"

#include <cctype>
#include <iostream>
#include <limits>

namespace fs = std::filesystem;

fasta_reader::fasta_reader(fs::path path, std::size_t buffer_size)
  : file{path}, path{path}, index{0} {
  if (!file.is_open()) {
    std::cerr << "Unable to open file, aborting...\n";
    exit(1);
  }
  buffer.resize(buffer_size - buffer_size%dna::size() + 1);
  load_buffer();
}

/**
 *  Advances the reader towards the next meaningful FASTA symbol.
 */
void fasta_reader::next_symbol() {
  index += dna::size();
  if (index >= buffer.size()-1) {
    if (!file.eof())
      load_buffer();
    else
      end_of_file = true;
  }
}

/**
 *  Loads the next data in the FASTA file into the current buffer.
 */
void fasta_reader::load_buffer() {
  auto position = 0lu;

  while (position < buffer.size()-1) {
    file.clear();
    if (file.peek() == '>' || file.peek() == '\n') {
      file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    file.getline(&buffer[position], buffer.size() - position);

    // When a newline is found: failbit == false, actual line size = gcount() - 1
    // When the buffer end is reached: failbit == true, actual line size = gcount()
    const auto size = file.fail() ? file.gcount() : file.gcount() - 1;
    position += size;
    
    if (file.eof()) {
      buffer.resize(position - position%dna::size());
      index = 0;
      return;
    }
  }
  index = 0;
}

/**
 * Merely an upper bound, as comments and newline characters are also included
 * in this count. As each byte character is a single base pair, the number of
 * base pairs is bounded by the size of the file in bytes.
 */
auto fasta_reader::size() -> std::size_t {
  return std::filesystem::file_size(path);
}

auto read_genome(const fs::path path) -> std::vector<dna> {
  if (!fs::is_regular_file(path)) {
    std::cerr << "Non-existent path, aborting...\n";
    exit(1);
  }

  auto result = std::vector<dna>{};
  auto file = fasta_reader{path};
  for (const auto&& element : file)
    result.emplace_back(element);

  return result;
}