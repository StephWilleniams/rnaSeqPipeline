#include <iostream>
#include <string>
#include <fstream>
#include <zlib.h>
#include <stdio.h>
#include <unistd.h> 
#include <chrono> 
#include "kseq.h"

// Initialize kseq for two separate files (R1 and R2)
KSEQ_INIT(gzFile, gzread)

// Helper function to write to a gzFile (Optimized with gzputs)
void write_fastq(gzFile out, kseq_t* seq, const std::string& umi = "", int slice_start = 0) {
    std::string header = seq->name.s;
    if (!umi.empty()) {
        header += "_" + umi; 
    }
    if (seq->comment.l > 0) {
        header += " " + std::string(seq->comment.s);
    }
    
    std::string out_seq = std::string(seq->seq.s).substr(slice_start);
    std::string out_qual = std::string(seq->qual.s).substr(slice_start);
    
    // Pre-build the string in memory, then blast it to the file using gzputs
    std::string record = "@" + header + "\n" + out_seq + "\n+\n" + out_qual + "\n";
    gzputs(out, record.c_str());
}

int main(int argc, char *argv[]) {
    std::string in_path1, in_path2, ext_path1, ext_path2, log_path, filtered_path1, filtered_path2;
    bool inputs_set = false, outputs_set = false;

    // Use getopt to parse standard flags
    int opt;
    while ((opt = getopt(argc, argv, "r:R:o:O:L:i:I:")) != -1) {
        switch (opt) {
            case 'r': in_path1 = optarg; break;
            case 'R': in_path2 = optarg; break;
            case 'o': ext_path1 = optarg; break;
            case 'O': ext_path2 = optarg; break;
            case 'L': log_path = optarg; break;
            case 'i': filtered_path1 = optarg; break;
            case 'I': filtered_path2 = optarg; break;
            default:
                std::cerr << "Usage: " << argv[0] << " -r <R1_in> -R <R2_in> [-o <ext_R1> -O <ext_R2>] [-L <log_file>] [-i <internal_R1> -I <internal_R2>]\n";
                return 1;
        }
    }

    // Validation: Ensure required inputs were provided
    if (in_path1.empty() || in_path2.empty()) {
        std::cerr << "Error: -r and -R are required.\n";
        return 1;
    }

    // Auto-generate outputs if not explicitly provided
    auto make_out_path = [](const std::string& path, const std::string& prefix) {
        size_t last_slash = path.find_last_of('/');
        return (last_slash == std::string::npos) ? prefix + path : path.substr(0, last_slash + 1) + prefix + path.substr(last_slash + 1);
    };

    if (ext_path1.empty()) ext_path1 = make_out_path(in_path1, "extracted_");
    if (ext_path2.empty()) ext_path2 = make_out_path(in_path2, "extracted_");
    
    // Auto-generate log path if not provided
    if (log_path.empty()) {
        log_path = make_out_path(in_path1, "extraction_progress_");
        size_t ext_pos = log_path.find(".fastq.gz");
        if (ext_pos != std::string::npos) log_path.replace(ext_pos, 9, ".log");
        else log_path += ".log";
    }

    // Always auto-generate internal paths based on the input names
    if (filtered_path1.empty()) filtered_path1 = make_out_path(in_path1, "internal_");
    if (filtered_path2.empty()) filtered_path2 = make_out_path(in_path2, "internal_");

    // Open actual processing files (Using "w1" for maximum compression speed)
    gzFile in_r1 = gzopen(in_path1.c_str(), "r");
    gzFile in_r2 = gzopen(in_path2.c_str(), "r");
    if (!in_r1 || !in_r2) {
        std::cerr << "Error: Could not open input files!\n";
        return 1;
    }
    kseq_t *seq1 = kseq_init(in_r1);
    kseq_t *seq2 = kseq_init(in_r2);

    gzFile ext_r1 = gzopen(ext_path1.c_str(), "w1");
    gzFile ext_r2 = gzopen(ext_path2.c_str(), "w1");
    gzFile int_r1 = gzopen(filtered_path1.c_str(), "w1");
    gzFile int_r2 = gzopen(filtered_path2.c_str(), "w1");
    
    // Open Log File
    std::ofstream log_file(log_path);

    const std::string anchor = "ATTGCGCAATG";
    int anchor_len = 11;
    long long read_count = 0;
    long long extracted_count = 0;
    const long long log_interval = 1000000; // Log every 1 million reads

    std::cout << "Processing Strict Regex UMI Extraction...\n";
    std::cout << "Outputs will be saved as:\n  " << ext_path1 << "\n  " << ext_path2 << "\n";
    std::cout << "Logging progress to:\n  " << log_path << "\n\n";

    log_file << "Starting Extraction Job...\n-----------------------------------\n";

    // --- START TIMER ---
    auto start_time = std::chrono::high_resolution_clock::now();

    while (kseq_read(seq1) >= 0 && kseq_read(seq2) >= 0) {
        read_count++;
        bool is_extracted = false;

        // --- PROGRESS LOGGER ---
        if (read_count % log_interval == 0) {
            log_file << "Progress: " << read_count << " reads processed...\n";
            log_file.flush(); 
            //std::cout << "Progress: " << read_count << " reads processed...\n"; 
        }

        int min_required_length = anchor_len + 8 + 3;

        // --- GREEDY BACKTRACKING BLOCK ---
        if (seq1->seq.l >= min_required_length) {
            int best_anchor_end = -1;
            int best_bio_start = -1;
            std::string best_umi = "";
            
            int search_limit = seq1->seq.l - min_required_length;
            
            // Sliding window to mimic {s<=2} (Strict Substitutions)
            for (int ii = 0; ii <= search_limit; ++ii) {
                int mismatches = 0;
                for (int jj = 0; jj < anchor_len; ++jj) {
                    if (seq1->seq.s[ii + jj] != anchor[jj]) {
                        mismatches++;
                    }
                }
                
                // If the anchor matches, IMMEDIATELY check the G requirement
                if (mismatches <= 2) {
                    int temp_anchor_end = ii + anchor_len - 1; 
                    int temp_bio_start = temp_anchor_end + 9;
                    
                    int g_count = 0;
                    while ((temp_bio_start + g_count) < seq1->seq.l && seq1->seq.s[temp_bio_start + g_count] == 'G') {
                        g_count++;
                    }

                    // EXACT REGEX BEHAVIOR: Must have >= 3 Gs to be a valid match
                    if (g_count >= 3) {
                        best_anchor_end = temp_anchor_end;
                        int discard_amount = (g_count > 5) ? 5 : g_count;
                        best_bio_start = temp_bio_start + discard_amount;
                        best_umi = std::string(seq1->seq.s).substr(best_anchor_end + 1, 8);
                        
                        // No break; statement, ensuring we grab the LAST valid match in the string
                    }
                }
            }

            // If we found at least one valid Anchor + G-stretch combination
            if (best_anchor_end != -1) {
                write_fastq(ext_r1, seq1, best_umi, best_bio_start);
                write_fastq(ext_r2, seq2, best_umi, 0);         
                
                extracted_count++;
                is_extracted = true;
            }
        }

        if (!is_extracted) {
            write_fastq(int_r1, seq1, "", 0);
            write_fastq(int_r2, seq2, "", 0);
        }
    }

    // --- STOP TIMER ---
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;

    // --- WRITE FINAL SUMMARY ---
    log_file << "-----------------------------------\n";
    log_file << "Job Complete!\n";
    log_file << "Total Pairs Evaluated: " << read_count << "\n";
    log_file << "Total UMIs Extracted:  " << extracted_count << "\n";
    log_file << "Time Elapsed:          " << elapsed_seconds.count() << " seconds\n";
    
    std::cout << "-----------------------------------\n";
    std::cout << "Job Complete!\n";
    std::cout << "Total Pairs Evaluated: " << read_count << "\n";
    std::cout << "Total UMIs Extracted:  " << extracted_count << "\n";
    std::cout << "Time Elapsed:          " << elapsed_seconds.count() << " seconds\n";

    kseq_destroy(seq1); kseq_destroy(seq2);
    gzclose(in_r1); gzclose(in_r2);
    gzclose(ext_r1); gzclose(ext_r2);
    gzclose(int_r1); gzclose(int_r2);
    log_file.close();

    return 0;
}