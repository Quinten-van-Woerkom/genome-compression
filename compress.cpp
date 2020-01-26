/**
 *  Entry point to the application.
 *  Compresses the given file using a canonicalized directed acyclic graph.
 */

#include <chrono>
#include <filesystem>
#include <iostream>
#include <vector>

#include "shared_tree.h"
#include "dna.h"
#include "fasta_reader.h"

void print_input(std::filesystem::path input_file, std::uintmax_t file_size) {
  std::cout
    << "\n============================================================\n"
    << " Input\n"
    << "============================================================\n"
    << " Filename:                  " << input_file << '\n'
    << " Size:                      " << bytes_to_string(file_size) << '\n'
    << " Nucleotides (upper bound): " << file_size << "\n\n";
}

void print_output(std::filesystem::path output_file, std::filesystem::path histogram,
  std::size_t compressed_size, std::size_t compressed_width, std::size_t file_size)
{

  std::cout
    << "\n============================================================\n"
    << " Output\n"
    << "============================================================\n";

  if (!output_file.empty())
    std::cout << " Filename:                  " << output_file << '\n';

  std::cout
    << " Size:                      " << bytes_to_string(compressed_size) << '\n'
    << " Nucleotides:               " << compressed_width*dna::size() << '\n'
    << " Compression ratio:         " << double(file_size)/double(compressed_size) << '\n';
  
  if (!histogram.empty())
    std::cout << " Histogram:                 " << histogram << '\n';
}

void print_tree_dimensions(shared_tree& tree, std::size_t width) {
  std::cout
    << "\n============================================================\n"
    << " Tree dimensions\n"
    << "============================================================\n"
    << " Leaf size:                 " << dna::size() << " nucleotides\n"
    << " Width:                     " << width << '\n'
    << " Depth:                     " << tree.depth() << '\n'
    << " Leaves:                    " << tree.leaf_count() << '\n'
    << " Nodes:                     " << tree.node_count() << '\n';
}

void print_timings(std::chrono::milliseconds construction, std::chrono::milliseconds sorting) {
  std::cout
    << "\n============================================================\n"
    << " Timings\n"
    << "============================================================\n"
    << " Tree construction:         " << construction.count() << " ms\n"
    << " Frequency sorting:         " << sorting.count() << " ms\n\n";
}

void print_statistics(std::size_t original_size, std::size_t compressed_size,
  std::size_t compressed_width, std::chrono::milliseconds construction,
  std::chrono::milliseconds sorting)
{
  std::cout << dna::size()
    << ',' << compressed_width
    << ',' << double(original_size)/double(compressed_size)
    << ',' << original_size
    << ',' << compressed_size
    << ',' << construction.count()
    << ',' << sorting.count()
    << ',' << construction.count() + sorting.count()
    << '\n';
}

void print_help() {
  std::cout
    << "Usage: compress [options] file...\n"
    << "Options:\n"
    << "\t--help\t\t\tPrints this documentation\n"
    << "\t--verbose\t\tPrint verbose output\n"
    << "\t--statistics\t\tPrint only numerical summary of output\n"
    << "\t--no-save\t\tDo not save the compressed file\n"
    << "\t--output=<file>\t\tWrite output to <file>, default being <input>.dag\n"
    << "\t--histogram=<file>\tSave histogram of node references in tree to <file>\n"
    << "\t--dna-size=<size>\tThe number of nucleotides stored per leaf node, default is 12\n";
}

auto parse_commands(int argc, char* argv[]) {
  std::filesystem::path input_file;
  std::filesystem::path output_file;
  std::filesystem::path histogram;
  bool verbose = false;
  bool statistics = false;
  bool save = true;
  std::size_t dna_size = 12;

  if (argc == 1) {
    std::cout << "Invalid command: argument <file> required.\n";
    std::cout << "Use --help for more information\n";
    exit(2);
  }

  for (auto i = 1u; i < argc; ++i) {
    auto argument = std::string_view{argv[i]};

    if (argument == "--help") {
      print_help();
      exit(0);
    } else if (argument == "--verbose") {
      verbose = true;
      continue;
    } else if (argument == "--statistics") {
      statistics = true;
      continue;
    } else if (argument.substr(0, 9) == "--output=") {
      argument.remove_prefix(9);
      output_file = argument;
      continue;
    } else if (argument.substr(0, 12) == "--histogram=") {
      argument.remove_prefix(12);
      histogram = argument;
      continue;
    } else if (argument == "--no-save") {
      save = false;
      continue;
    } else if (argument.substr(0, 11) == "--dna-size=") {
      argument.remove_prefix(11);
      std::cout << argument << '\n';
      dna_size = std::atoi(argument.data());
      continue;
    } else { // Interpret as name of input file
      if (!input_file.empty()) {
        std::cout << "Compression of multiple files at once is currently not supported.\n";
        exit(1);
      }
      input_file = argument;
    }
  }

  if (verbose && statistics) {
    std::cout << "Invalid flag combination: --verbose and --statistics are mutually exclusive\n";
    std::cout << "Use --help for more information\n";
    exit(2);
  }

  if (input_file.empty()) {
    std::cout << "Invalid command: argument <file> required.\n";
    std::cout << "Use --help for more information\n";
    exit(2);
  }

  if (output_file.empty() && save) {
    output_file = input_file;
    output_file.replace_extension(".dag");
  }

  return std::tuple{input_file, output_file, histogram, verbose, statistics, dna_size};
}

int main(int argc, char* argv[]) {
  auto [input_file, output_file, histogram, verbose, statistics, dna_size] = parse_commands(argc, argv);
  dna::size(dna_size);

  if (!std::filesystem::is_regular_file(input_file)) {
    std::cout << "Invalid filename: " << input_file << '\n';
    exit(2);
  }

  auto original_size = std::filesystem::file_size(input_file);

  if (verbose)
    print_input(input_file, original_size);


  auto start = std::chrono::high_resolution_clock::now();
  auto compressed = shared_tree{input_file, verbose};
  auto end = std::chrono::high_resolution_clock::now();
  auto construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  
  start = std::chrono::high_resolution_clock::now();
  compressed.sort_tree(verbose);
  end = std::chrono::high_resolution_clock::now();
  auto sorting_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

  auto compressed_size = compressed.bytes();
  auto compressed_width = compressed.width();


  if (!histogram.empty())
    compressed.store_histogram(histogram);

  if (!output_file.empty())
    compressed.save(output_file);

  if (verbose) {
    print_output(output_file, histogram, compressed_size, compressed_width, original_size);
    print_tree_dimensions(compressed, compressed_width);
    print_timings(construction_time, sorting_time);
  }

  if (statistics) {
    print_statistics(original_size, compressed_size, compressed_width, construction_time, sorting_time);
  }

  return 0;
}