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
  : file{path}, path{path} {
  if (!file.is_open()) {
    std::cerr << "Unable to open file, aborting...\n";
    exit(1);
  }


  // Make sure that we do not allocate an unnecessarily big buffer.
  const auto file_size = std::filesystem::file_size(path);
  const auto file_strands = file_size/dna::size()+dna::size();

  if (file_strands < buffer_size) {
    buffer.resize(file_strands);
    char_buffer.resize(file_strands*dna::size() + 1);
  } else {
    buffer.resize(buffer_size);
    char_buffer.resize(buffer_size*dna::size() + 1);
  }
  background_loader = std::thread{&fasta_reader::load_buffer, this};
}

/**
 *  Loads the next data in the FASTA file into the background buffer.
 */
void fasta_reader::load_buffer() {
  if (file.eof()) {
    end_of_file = true;
    return;
  }

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
      char_buffer.resize((position+1)/dna::size() * dna::size());
      buffer.resize(char_buffer.size() / dna::size());
      break;
    }
  }

  for (auto i = 0u; i < buffer.size(); ++i)
    buffer[i] = dna{std::string_view{&char_buffer[i*dna::size()], dna::size()}};
}

/**
 * Merely an upper bound, as comments and newline characters are also included
 * in this count. As each byte character is a single base pair, the number of
 * base pairs is bounded by the size of the file in bytes.
 */
auto fasta_reader::size() -> std::size_t {
  return std::filesystem::file_size(path);
}

/**
 * Reads the data currently in the buffer into the vector passed.
 * Does this through a move swap, so that the data itself is not actually read.
 * Returns true if the read was successful, false otherwise.
 */
bool fasta_reader::read_into(std::vector<dna>& vector) {
  if (end_of_file) return false;
  if (background_loader.joinable()) background_loader.join();
  std::swap(buffer, vector);

  if (!file.eof()) {
    buffer.resize(vector.size());
    background_loader = std::thread{&fasta_reader::load_buffer, this};
  } else {
    end_of_file = true;
  }

  return vector.size() != 0;
}

auto read_genome(const std::filesystem::path path) -> std::vector<dna> {
  if (!std::filesystem::is_regular_file(path)) {
    std::cerr << "Non-existent path, aborting...\n";
    exit(1);
  }

  auto result = std::vector<dna>{};
  auto buffer = std::vector<dna>{};
  auto file = fasta_reader{path};
  while (file.read_into(buffer)) {
    for (const auto& element : buffer)
      result.emplace_back(element);
  }

  return result;
}