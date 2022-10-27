#pragma once
#include <iostream>
#include <set>
#include <string>

extern "C" {
#include "taf.h"
}

int uce_main(int argc, char** argv);

void compute_taf_uces(FILE* input, std::ostream& os, int64_t min_len, int64_t min_depth, const std::set<std::string>& eclusion_set);

// divide by first occurrence of separator
std::pair<std::string, std::string> parse_sample_contig(const char* sequence_name, char separator = '.');

// get the coverage at a certain column
int64_t column_depth(Alignment* alignment, int64_t column);
