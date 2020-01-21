/**
 *  Buffered file reader implementation for single-FASTA DNA sequences.
 *  Allows one to read bigger DNA sequences without requiring them to be fully
 *  loaded in memory.
 */

#include "fasta_reader.h"
#include "dna.h"

#include <iostream>
#include <limits>

fasta_reader::fasta_reader(std::filesystem::path path, std::size_t buffer_size)
  : file{path}, path{path}, index{0} {
  if (!file.is_open()) {
    std::cerr << "Unable to open file, aborting...\n";
    exit(1);
  }

  // Make sure that we do not allocate an unnecessarily big buffer.
  const auto file_size = std::filesystem::file_size(path);
  if (file_size < buffer_size*(dna::size()/2))
    buffer_size = file_size;

  buffer.resize(buffer_size);
  background_buffer.resize(buffer_size);
  char_buffer.resize(buffer_size*dna::size() + 1);
  load_buffer();
  std::swap(background_buffer, buffer);
  background_loader = std::thread{&fasta_reader::load_buffer, this};
}

/**
 *  Advances the reader towards the next meaningful FASTA symbol.
 */
void fasta_reader::next_symbol() {
  ++index;
  if (index >= buffer.size()) {
    if (!file.eof()) swap_buffers();
    else end_of_file = true;
  }
}

/**
 * Waits until the background thread has finished loading the next buffer, and
 * then swaps it with the current buffer, afterwards dispatching the background
 * thread to load the next buffer.
 */
void fasta_reader::swap_buffers() {
  if (background_loader.joinable()) background_loader.join();
  std::swap(background_buffer, buffer);
  index = 0;
  background_loader = std::thread{&fasta_reader::load_buffer, this};
}

/**
 *  Loads the next data in the FASTA file into the background buffer.
 */
void fasta_reader::load_buffer() {
  auto position = 0lu;
  while (position < char_buffer.size()-1) {
    file.clear();
    if (file.peek() == '>' || file.peek() == '\n') {
      file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    file.getline(&char_buffer[position], char_buffer.size() - position);

    // When a newline is found: failbit == false, actual line size = gcount() - 1
    // When the buffer end is reached: failbit == true, actual line size = gcount()
    const auto size = file.fail() ? file.gcount() : file.gcount() - 1;
    position += size;
    
    if (file.eof()) {
      char_buffer.resize(position+1 - (position+1)%dna::size());
      background_buffer.resize(char_buffer.size() / dna::size());
      break;
    }
  }

  for (auto i = 0u; i < background_buffer.size(); ++i)
    background_buffer[i] = dna{std::string_view{&char_buffer[i*dna::size()], dna::size()}};
}

/**
 * Merely an upper bound, as comments and newline characters are also included
 * in this count. As each byte character is a single base pair, the number of
 * base pairs is bounded by the size of the file in bytes.
 */
auto fasta_reader::size() -> std::size_t {
  return std::filesystem::file_size(path);
}

auto read_genome(const std::filesystem::path path) -> std::vector<dna> {
  if (!std::filesystem::is_regular_file(path)) {
    std::cerr << "Non-existent path, aborting...\n";
    exit(1);
  }

  auto result = std::vector<dna>{};
  auto file = fasta_reader{path};
  for (const auto&& element : file)
    result.emplace_back(element);

  return result;
}