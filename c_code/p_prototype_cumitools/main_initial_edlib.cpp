#include <iostream>
#include <string>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "edlib.h"

// Initialize kseq for two separate files (R1 and R2)
KSEQ_INIT(gzFile, gzread)

// Helper function to write to a gzFile
void write_fastq(gzFile out, kseq_t* seq, const std::string& umi = "", int slice_start = 0) {
    // 1. Grab the first part of the header (before the space)
    std::string header = seq->name.s;
    
    // 2. Append the UMI if we extracted one
    if (!umi.empty()) {
        header += "_" + umi; 
    }
    
    // 3. CRITICAL FIX: If the header had a second part (after the space), paste it back!
    if (seq->comment.l > 0) {
        header += " " + std::string(seq->comment.s);
    }
    
    // Slice sequence and quality strings
    std::string out_seq = std::string(seq->seq.s).substr(slice_start);
    std::string out_qual = std::string(seq->qual.s).substr(slice_start);
    
    // Write out the rebuilt read
    gzprintf(out, "@%s\n%s\n+\n%s\n", header.c_str(), out_seq.c_str(), out_qual.c_str());
}

int main(int argc, char *argv[]) {
    // Now we only require 3 arguments: Program Name, R1_in, R2_in
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <R1_in.fastq.gz> <R2_in.fastq.gz>\n";
        return 1;
    }

    std::string in_path1 = argv[1];
    std::string in_path2 = argv[2];

    // Lambda function to safely insert a prefix into the filename, ignoring directories
    auto make_out_path = [](const std::string& path, const std::string& prefix) {
        size_t last_slash = path.find_last_of('/');
        if (last_slash == std::string::npos) {
            return prefix + path; // No directories, just prepend
        } else {
            // Keep the directory, but prepend the prefix to the actual filename
            return path.substr(0, last_slash + 1) + prefix + path.substr(last_slash + 1);
        }
    };

    // Automatically generate the 4 output filenames
    std::string ext_path1 = make_out_path(in_path1, "extracted_");
    std::string ext_path2 = make_out_path(in_path2, "extracted_");
    std::string int_path1 = make_out_path(in_path1, "internal_");
    std::string int_path2 = make_out_path(in_path2, "internal_");

    // Open Inputs
    gzFile in_r1 = gzopen(in_path1.c_str(), "r");
    gzFile in_r2 = gzopen(in_path2.c_str(), "r");
    if (!in_r1 || !in_r2) {
        std::cerr << "Error: Could not open input files!\n";
        return 1;
    }
    kseq_t *seq1 = kseq_init(in_r1);
    kseq_t *seq2 = kseq_init(in_r2);

    // Open Outputs using our generated filenames
    gzFile ext_r1 = gzopen(ext_path1.c_str(), "w");
    gzFile ext_r2 = gzopen(ext_path2.c_str(), "w");
    gzFile int_r1 = gzopen(int_path1.c_str(), "w");
    gzFile int_r2 = gzopen(int_path2.c_str(), "w");

    const char* anchor = "ATTGCGCAATG";
    int anchor_len = 11;
    int read_count = 0;
    int extracted_count = 0;

    std::cout << "Processing UMI Extraction...\n";
    std::cout << "Outputs will be saved as:\n  " << ext_path1 << "\n  " << ext_path2 << "\n\n";

    // Loop through R1 and R2 synchronously
    while (kseq_read(seq1) >= 0 && kseq_read(seq2) >= 0) {
        read_count++;
        bool is_extracted = false;

        // Ensure R1 is long enough to even contain an anchor + UMI
        if (seq1->seq.l >= 30) {
            
            // Search the first 30 bases of R1 for the anchor
            std::string search_window = std::string(seq1->seq.s).substr(0, 30);
            
            // HW mode: finds sequence within the target, allowing 2 mismatches/indels
            EdlibAlignResult result = edlibAlign(anchor, anchor_len, search_window.c_str(), search_window.length(),
                                                 edlibNewAlignConfig(2, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0));

            if (result.status == EDLIB_STATUS_OK && result.numLocations > 0) {
                int anchor_end = result.endLocations[0]; // 0-based index of where anchor finishes

                // Check bounds: ensure there's enough room for an 8-base UMI
                if (anchor_end + 8 < seq1->seq.l) {
                    std::string umi = std::string(seq1->seq.s).substr(anchor_end + 1, 8);
                    
                    int bio_start = anchor_end + 9;
                    
                    // Conditionally discard the Poly-G stretch (biological shift)
                    // Mimics regex G{3,5}: Requires at least 3 Gs, discards max of 5.
                    int g_count = 0;
                    while ((bio_start + g_count) < seq1->seq.l && seq1->seq.s[bio_start + g_count] == 'G') {
                        g_count++;
                    }

                    if (g_count >= 3) {
                        // Discard the Gs (up to a maximum of 5)
                        int discard_amount = (g_count > 5) ? 5 : g_count;
                        bio_start += discard_amount;
                    }
                    // If g_count is 0, 1, or 2, bio_start remains completely unchanged.

                    // Write to Extracted Files with modified headers
                    write_fastq(ext_r1, seq1, umi, bio_start); // R1 gets sliced
                    write_fastq(ext_r2, seq2, umi, 0);         // R2 gets untouched sequence
                    
                    extracted_count++;
                    is_extracted = true;
                }
            }
            edlibFreeAlignResult(result);
        }

        // If it failed the search, write to Internal files unmodified
        if (!is_extracted) {
            write_fastq(int_r1, seq1, "", 0);
            write_fastq(int_r2, seq2, "", 0);
        }
    }

    std::cout << "Complete! Evaluated: " << read_count << " pairs. Extracted: " << extracted_count << " UMIs.\n";

    // Cleanup
    kseq_destroy(seq1); kseq_destroy(seq2);
    gzclose(in_r1); gzclose(in_r2);
    gzclose(ext_r1); gzclose(ext_r2);
    gzclose(int_r1); gzclose(int_r2);

    return 0;
}