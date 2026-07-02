import itertools

# Your provided list of headers
raw_headers = [
    "1  - @VH01192:275:AAJCGLVM5:1:1101:21678:1000 1:N:0:CTGCGGTT+CTAGCGCT",
    "2  - @VH01192:275:AAJCGLVM5:1:1101:19027:1000 1:N:0:TTGTAACC+TCGATATT",
    "3  - @VH01192:275:AAJCGLVM5:1:1101:19064:1000 1:N:0:GGGCTTGG+CGTCTGCG",
    "4  - @VH01192:275:AAJCGLVM5:1:1101:19083:1019 1:N:0:AGGTCCAA+TACTCATA",
    "5  - @VH01192:275:AAJCGLVM5:1:1101:21772:1019 1:N:0:ATGCACTG+ACGCACCT",
    "6  - @VH01192:275:AAJCGLVM5:1:1101:19708:1000 1:N:0:GTTTGTCA+GTATGTTT",
    "7  - @VH01192:275:AAJCGLVM5:1:1101:18761:1038 1:N:0:CAAGCTAG+CGCTATGT",
    "8  - @VH01192:275:AAJCGLVM5:1:1101:22435:1000 1:N:0:TGGATTGA+TATCGCAT",
    "9  - @VH01192:275:AAJCGLVM5:1:1101:18913:1000 1:N:0:AGTCCAGG+TCTGTTAG",
    "10 - @VH01192:275:AAJCGLVM5:1:1101:18534:1038 1:N:0:GACCTGAA+CTCACCAA",
    "11 - @VH01192:275:AAJCGLVM5:1:1101:18591:1019 1:N:0:TTTCTACT+GAACCGCG",
    "12 - @VH01192:275:AAJCGLVM5:1:1101:19424:1019 1:N:0:CTTTCGTC+AGGTTATA"
]

def get_hamming_distance(seq1, seq2):
    """Calculates the number of differing characters between two strings."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length.")
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def main():
    # 1. Extract just the barcodes and map them to their sample number
    barcodes = {}
    for line in raw_headers:
        # Split by spaces to isolate the Sample ID and the actual header string
        parts = line.strip().split()
        sample_id = parts[0]
        header_data = parts[-1]
        
        # Split the header by colons to grab the very last section (the barcode)
        barcode = header_data.split(':')[-1]
        barcodes[sample_id] = barcode

    # 2. Compare all unique pairs
    min_distance = float('inf')
    closest_pairs = []

    # itertools.combinations ensures we don't compare a sample to itself 
    # and doesn't repeat comparisons (e.g., checks 1 vs 2, but skips 2 vs 1)
    for (id1, bc1), (id2, bc2) in itertools.combinations(barcodes.items(), 2):
        dist = get_hamming_distance(bc1, bc2)
        
        if dist < min_distance:
            min_distance = dist
            closest_pairs = [(id1, id2, dist)]
        elif dist == min_distance:
            closest_pairs.append((id1, id2, dist))

    # 3. Print the summary report
    print("=== Barcode Edit Distance Report ===")
    print(f"Total Barcodes Analysed: {len(barcodes)}")
    print(f"Minimum Edit Distance Found: {min_distance}\n")
    
    print("Closest Sample Pairs (Highest Risk of Misassignment):")
    for id1, id2, dist in closest_pairs:
        print(f"  Sample {id1} vs Sample {id2} -> {dist} mismatches")
        print(f"    {id1}: {barcodes[id1]}")
        print(f"    {id2}: {barcodes[id2]}\n")

    # 4. Assess the safety of a 2-mismatch rescue strategy
    print("=== Strategy Assessment ===")
    if min_distance > 4:
        print("Status: SAFE.")
        print("Allowing up to 2 substitutions per read is mathematically safe.")
        print("Even with 2 errors, a read cannot accidentally match another sample's barcode.")
    elif min_distance == 4:
        print("Status: BORDERLINE.")
        print("Allowing up to 2 substitutions could result in a tie (a read might be exactly 2 mismatches away from TWO different samples).")
        print("You must implement logic to discard ties.")
    else:
        print("Status: UNSAFE.")
        print(f"Because the minimum distance is {min_distance}, allowing 2 mismatches will cause cross-contamination between samples.")

if __name__ == "__main__":
    main()