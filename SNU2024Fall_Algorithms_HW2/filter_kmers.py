import math
import sys

def read_kmer():
    for line in sys.stdin:
        kmer, count = line.strip().split(',')
        yield (kmer, int(count))
        
def get_cutoff_multiplier(k):
    if k > 20:
        return 0.01
    elif 13 <= k <= 20:
        return 0.10
    elif k == 12:
        return 0.20
    elif k == 11:
        return 0.70
    elif k == 10:
        return 0.80
    elif k in [8, 9]:
        return 0.90
    elif k in [6, 7]:
        return 0.96
    elif 3 <= k <= 5:
        return 0.95
    elif k == 2:
        return 0.90
    elif k == 1:
        return 0.50
    else:
        # Default multiplier if k doesn't match any condition
        return 1.0

def filter_kmers(kmer_generator, k):
    filtered_kmers = []
    prefix_set = set()
    suffix_set = set()
    
    min_overlap = math.ceil(k / 2)
    
    kmer_list = list(kmer_generator)
    cutoff = int(len(kmer_list) * get_cutoff_multiplier(k))
    kmer_list = kmer_list[cutoff:]
    
    for kmer, _ in kmer_list:
        overlap = False
        for overlap_len in range(min_overlap, k):
            prefix = kmer[:overlap_len]
            suffix = kmer[-overlap_len:]
            
            if suffix in prefix_set or prefix in suffix_set:
                overlap = True
                break
        
        if not overlap:
            filtered_kmers.append(kmer)       
            prefixes = [kmer[:ol] for ol in range(min_overlap, k)]
            suffixes = [kmer[-ol:] for ol in range(min_overlap, k)]
            prefix_set.update(prefixes)
            suffix_set.update(suffixes)
    
    return filtered_kmers

def main(k):
    kmer_generator = read_kmer()
    filtered_kmers = filter_kmers(kmer_generator, k)
    for kmer in filtered_kmers:
        print(kmer)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python filter_kmers.py <k>")
        sys.exit(1)
    k = int(sys.argv[1])
    main(k)