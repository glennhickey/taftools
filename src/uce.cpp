#include "uce.hpp"
#include <unistd.h>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <cassert>
#include <algorithm>
#include <vector>
#include <sys/stat.h>
#include <set>

using namespace std;


static void help(char** argv) {
  cerr << "usage: taftools " << argv[0] << " [options]" << endl
       << "Compute Ultra-Conserved Elements from a TAF file" << endl
       << "input options: " << endl
       << "    -i, --input-file FILE                   Input file in TAF format [stdin]" << endl
       << "    -m, --min-length N                      Output regions must be at least N bp long [required]" << endl
       << "    -d, --min-depth F                       Output regions must be conserved in N samples [required]" << endl
       << "    -x, --exclude-sample                    Columns containing given sample are not considered conserved (multiple allowed)" << endl
       << endl;
}    

int uce_main(int argc, char** argv) {

    string input_taf_path;
    int64_t min_uce_len = 0;
    int64_t min_uce_depth = 0;
    set<string> exclusion_set;
    
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"input-file", required_argument, 0, 'i'},
            {"min-length", required_argument, 0, 'm'},
            {"min-depth", required_argument, 0, 'd'},
            {"exclude-sample", required_argument, 0, 'x'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hi:m:d:x:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'i':
            input_taf_path = optarg;
            break;
        case 'm':
            min_uce_len = stol(optarg);
            break;
        case 'd':
            min_uce_depth = stol(optarg);
            break;
        case 'x':
            exclusion_set.insert(optarg);
            break;
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (min_uce_depth <= 0) {
        cerr << "[uce] Depth threshold must be specified with -d" << endl;
        help(argv);
        return 1;
    }
    if (min_uce_len <= 0) {
        cerr << "[uce] Length threshold must be specified with -l" << endl;
        help(argv);
        return 1;
    }

    FILE *input = input_taf_path.empty() ? stdin : fopen(input_taf_path.c_str(), "r");
    if (!input) {
        cerr << "[uce] Unable to open input file " << input_taf_path << endl;
        return 1;
    }

    // compute the uces, printing to stdout
    compute_taf_uces(input, cout, min_uce_len, min_uce_depth, exclusion_set);
    
    if (input != NULL) {
        fclose(input);
    }
    
    return 0;
}

void compute_taf_uces(FILE* input, std::ostream& os, int64_t min_len, int64_t min_depth, const set<string>& exclusion_set) {
    // make the iterator
    LI *li = LI_construct(input);

    // read the TAF header
    stList *tags = taf_read_header(li);
    assert(stList_length(tags) % 2 == 0);

    // keep track of curent uce
    string uce_contig;
    int64_t uce_start = -1;
    int64_t uce_last = -1;
    int64_t uce_mindepth = 0;

    // keep track of sample names
    vector<string> samples;
    
    // scan the blocks
    bool run_length_encode_bases = false;
    Alignment *alignment, *p_alignment = NULL;
    while((alignment = taf_read_block(p_alignment, run_length_encode_bases, li)) != NULL) {

        bool skip = false;
        // don't waste time
        if (alignment->row_number < min_depth) {
            skip = true;
        }

        if (p_alignment && p_alignment->row_number > 0 &&
            strcmp(alignment->row->sequence_name, p_alignment->row->sequence_name) != 0) {
            // should solidify this with an explicit reference name check
            skip = true;
        }

        if (alignment->row->strand == 0) {
            cerr << "[uce] negative reference strand found!" << endl;
            exit(1);
        }

        if (!skip) { 
            // make a sample set (warning this is getting pretty inefficient)
            unordered_set<string> sample_set;
        
            // parse the samples
            samples.clear();
            for (Alignment_Row* row = alignment->row; row != NULL; row = row->n_row) {
                samples.push_back(parse_sample_contig(row->sequence_name).first);
                sample_set.insert(samples.back());
            }

            // now iterate the columns
            int64_t ref_offset = alignment->row->start;
            for (int64_t col = 0; col < alignment->row->length; ++col) {
                if (alignment->row->bases[col] != '-') {                
                    int64_t column_depth = 0;            
                    if (sample_set.size() == samples.size() && exclusion_set.empty()) {
                        // if we have no duplications, it's a simple count
                        for (Alignment_Row* row = alignment->row->n_row; row != NULL; row = row->n_row) {
                            if (toupper(row->bases[col]) == toupper(alignment->row->bases[col])) {
                                ++column_depth;
                            }
                        }
                    } else {
                        // we have duplications, then we do a much slower set comparison
                        // todo: something more clever (an easy thing would be to track which rows are potential dupes)
                        sample_set.clear();
                        sample_set.insert(samples[0]);
                        int64_t row_i = 1;
                        for (Alignment_Row* row = alignment->row->n_row; row != NULL; row = row->n_row, ++row_i) {
                            if (toupper(row->bases[col]) == toupper(alignment->row->bases[col])) {
                                string& sample = samples[row_i];
                                if (exclusion_set.count(sample)) {
                                    // we found an excluded sample, the column is void
                                    column_depth = 0;
                                    break;
                                }
                                if (!sample_set.count(sample)) {
                                    ++column_depth;
                                    sample_set.insert(sample);
                                }
                            }
                        }
                    }
                    bool continue_uce = column_depth >= min_depth && ref_offset == uce_last + 1 && uce_start >= 0 &&
                        uce_contig == alignment->row->sequence_name;

                    if (continue_uce) {
                        uce_last = ref_offset;
                        uce_mindepth = std::min(uce_mindepth, column_depth);
                    } else {
                        // gotta close the current uce
                        if (uce_last - uce_start >= min_len - 1) {
                            os << uce_contig << "\t" << uce_start << "\t" << (uce_last + 1) << "\t" << uce_mindepth << "\n";
                        }
                        if (column_depth >= min_depth) {
                            // we start a new one.
                            uce_contig = alignment->row->sequence_name;
                            uce_start = ref_offset;
                            uce_last = ref_offset;
                            uce_mindepth = column_depth;
                        } else {
                            // we have no open uce
                            uce_start = -1;
                            uce_last = -1;
                        }
                    }
                    ++ref_offset;
                }
            }
        }
                        
        // Clean up the previous alignment
        if(p_alignment != NULL) {
            alignment_destruct(p_alignment);
        }
        p_alignment = alignment; // Update the previous alignment
    }
    if(p_alignment != NULL) { // Clean up the final alignment
        alignment_destruct(p_alignment);
    }

    // there could be an open uce
    if (uce_last - uce_start >= min_len - 1) {
        os << uce_contig << "\t" << uce_start << "\t" << (uce_last + 1) << "\t" << uce_mindepth << "\n";
    }                        
    
    // clean up
    LI_destruct(li);    
}

pair<string, string> parse_sample_contig(const char* sequence_name, char separator) {
    const char* match = strchr(sequence_name, separator);
    if (match == NULL) {
        return make_pair("NULL", sequence_name);
    }
    pair<string, string> sc;

    size_t slen = match - sequence_name;
    sc.first.reserve(slen);
    for (size_t i = 0; i < slen; ++i) {
        sc.first.push_back(sequence_name[i]);
    }
    sc.second = match + 1;
    return sc;
}
