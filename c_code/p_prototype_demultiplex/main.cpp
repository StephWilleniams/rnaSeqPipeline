#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <zlib.h>
#include <stdio.h>
#include <unistd.h> 
#include <chrono> 
#include <iomanip>
#include "kseq.h"

// Initialise kseq for two separate files (R1 and R2)
KSEQ_INIT(gzFile, gzread)

// Helper function to calculate Hamming distance between two strings
int calculate_mismatches(const std::string& read_barcode, const std::string& known_barcode) {
    if (read_barcode.length() != known_barcode.length()) return 999;
    
    int mismatches = 0;
    for (size_t i = 0; i < read_barcode.length(); ++i) {
        if (read_barcode[i] != known_barcode[i]) {
            mismatches++;
        }
    }
    return mismatches;
}

// Helper function to write to a gzFile (Optimised with gzputs)
void write_fastq(gzFile out, kseq_t* seq) {
    std::string header = seq->name.s;
    if (seq->comment.l > 0) {
        header += " " + std::string(seq->comment.s);
    }
    
    // Pre-build the string in memory, then blast it to the file using gzputs
    std::string record = "@" + header + "\n" + std::string(seq->seq.s) + "\n+\n" + std::string(seq->qual.s) + "\n";
    gzputs(out, record.c_str());
}

int main(int argc, char *argv[]) {
    std::string in_path1, in_path2, out_prefix, log_path;
    int max_mismatches = 3; // Default to a tolerance of 3

    // Parse standard flags
    int opt;
    while ((opt = getopt(argc, argv, "r:R:o:L:m:")) != -1) {
        switch (opt) {
            case 'r': in_path1 = optarg; break;
            case 'R': in_path2 = optarg; break;
            case 'o': out_prefix = optarg; break;
            case 'L': log_path = optarg; break;
            case 'm': max_mismatches = std::stoi(optarg); break;
            default:
                std::cerr << "Usage: " << argv[0] << " -r <R1_in> -R <R2_in> -o <output_prefix> [-m <max_mismatches>] [-L <log_file>]\n";
                return 1;
        }
    }

    if (in_path1.empty() || in_path2.empty() || out_prefix.empty()) {
        std::cerr << "Error: -r, -R, and -o (output prefix) are required.\n";
        return 1;
    }

    if (log_path.empty()) {
        log_path = out_prefix + "demultiplex_progress.log";
    }

    // Define the 12 known barcodes
    const std::vector<std::string> known_barcodes = {
        "CTGCGGTT+CTAGCGCT", "TTGTAACC+TCGATATT", "GGGCTTGG+CGTCTGCG",
        "AGGTCCAA+TACTCATA", "ATGCACTG+ACGCACCT", "GTTTGTCA+GTATGTTT",
        "CAAGCTAG+CGCTATGT", "TGGATTGA+TATCGCAT", "AGTCCAGG+TCTGTTAG",
        "GACCTGAA+CTCACCAA", "TTTCTACT+GAACCGCG", "CTTTCGTC+AGGTTATA"
    };
    const int num_samples = known_barcodes.size();

    // Open input streams
    gzFile in_r1 = gzopen(in_path1.c_str(), "r");
    gzFile in_r2 = gzopen(in_path2.c_str(), "r");
    if (!in_r1 || !in_r2) {
        std::cerr << "Error: Could not open input files!\n";
        return 1;
    }
    kseq_t *seq1 = kseq_init(in_r1);
    kseq_t *seq2 = kseq_init(in_r2);

    // Initialise arrays for output streams and tracking
    std::vector<gzFile> out_r1(num_samples + 1);
    std::vector<gzFile> out_r2(num_samples + 1);
    std::vector<long long> sample_counts(num_samples + 1, 0); // Index 12 is Undetermined

    // Open output streams (12 samples + 1 undetermined)
    for (int i = 0; i < num_samples; ++i) {
        std::string f1 = out_prefix + "sample_" + std::to_string(i + 1) + "_R1.fastq.gz";
        std::string f2 = out_prefix + "sample_" + std::to_string(i + 1) + "_R2.fastq.gz";
        out_r1[i] = gzopen(f1.c_str(), "w1");
        out_r2[i] = gzopen(f2.c_str(), "w1");
    }
    std::string un_f1 = out_prefix + "undetermined_R1.fastq.gz";
    std::string un_f2 = out_prefix + "undetermined_R2.fastq.gz";
    out_r1[num_samples] = gzopen(un_f1.c_str(), "w1");
    out_r2[num_samples] = gzopen(un_f2.c_str(), "w1");

    std::ofstream log_file(log_path);
    long long read_count = 0;
    long long failed_qc_count = 0;
    const long long log_interval = 1000000;

    std::cout << "Processing Custom Demultiplexing...\n";
    std::cout << "Max Mismatches Tolerated: " << max_mismatches << "\n";
    std::cout << "Outputs will be prefixed with: " << out_prefix << "\n";
    std::cout << "Logging progress to: " << log_path << "\n\n";

    log_file << "Starting Demultiplexing Job...\n-----------------------------------\n";

    auto start_time = std::chrono::high_resolution_clock::now();

    while (kseq_read(seq1) >= 0 && kseq_read(seq2) >= 0) {
        read_count++;

        if (read_count % log_interval == 0) {
            log_file << "Progress: " << read_count << " read pairs processed...\n";
            log_file.flush(); 
        }

        std::string comment = seq1->comment.s ? seq1->comment.s : "";
        
        // 1. Check vendor QC Flag
        if (comment.length() > 2 && comment[2] == 'Y') {
            failed_qc_count++;
            write_fastq(out_r1[num_samples], seq1);
            write_fastq(out_r2[num_samples], seq2);
            sample_counts[num_samples]++;
            continue;
        }

        // 2. Extract barcode from the end of the comment
        size_t last_colon = comment.find_last_of(':');
        std::string read_barcode = "";
        if (last_colon != std::string::npos && last_colon + 1 < comment.length()) {
            read_barcode = comment.substr(last_colon + 1);
        }

        int best_match_idx = -1;
        int min_dist = 999;
        int ties = 0;

        // 3. Calculate Hamming distances against all 12 known barcodes
        if (read_barcode.length() == 17) { // 8 bases + '+' + 8 bases
            for (int i = 0; i < num_samples; ++i) {
                int dist = calculate_mismatches(read_barcode, known_barcodes[i]);
                if (dist < min_dist) {
                    min_dist = dist;
                    best_match_idx = i;
                    ties = 1;
                } else if (dist == min_dist) {
                    ties++;
                }
            }
        }

        // 4. Assign read based on tolerance and lack of ties
        int target_idx = num_samples; // Default to undetermined
        if (best_match_idx != -1 && min_dist <= max_mismatches && ties == 1) {
            target_idx = best_match_idx;
        }

        // 5. Write to appropriate stream
        write_fastq(out_r1[target_idx], seq1);
        write_fastq(out_r2[target_idx], seq2);
        sample_counts[target_idx]++;
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;

    // --- WRITE FINAL SUMMARY ---
    auto write_summary = [&](std::ostream& os) {
        os << "-----------------------------------\n";
        os << "Job Complete!\n";
        os << "Time Elapsed:          " << elapsed_seconds.count() << " seconds\n";
        os << "Total Pairs Evaluated: " << read_count << "\n";
        os << "Vendor QC Failures:    " << failed_qc_count << "\n\n";
        
        os << "Breakdown by Sample:\n";
        long long total_assigned = 0;
        for (int i = 0; i < num_samples; ++i) {
            os << "  Sample " << std::setw(2) << (i + 1) << ": " << sample_counts[i] << " reads\n";
            total_assigned += sample_counts[i];
        }
        os << "\nTotal Assigned:        " << total_assigned << " (" << std::fixed << std::setprecision(2) << (double)total_assigned/read_count * 100 << "%)\n";
        os << "Total Undetermined:    " << sample_counts[num_samples] << "\n";
    };

    write_summary(log_file);
    write_summary(std::cout);

    // Clean up
    kseq_destroy(seq1); kseq_destroy(seq2);
    gzclose(in_r1); gzclose(in_r2);
    for (int i = 0; i <= num_samples; ++i) {
        gzclose(out_r1[i]);
        gzclose(out_r2[i]);
    }
    log_file.close();

    return 0;
}