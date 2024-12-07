import sys
from collections import Counter
import csv
from collections import defaultdict 

encode_array = [-1] * 128
encode_array[65] = 0b00  
encode_array[67] = 0b01  
encode_array[71] = 0b10  
encode_array[84] = 0b11 

decode_list = ['A', 'C', 'G', 'T']

def main():
    if len(sys.argv) != 3:
        print("Usage: python kmer_counter_final.py [k-mer] [input_file]")
        sys.exit(1)
    k = int(sys.argv[1])
    input_file = sys.argv[2].strip()
    student_id = '202116621'
    output_file = f"{student_id}.txt"
    kmer_counts = process_file(k, input_file)
    write_output(kmer_counts, output_file, k)

def process_file(k, input_file):
    kmer_counts = Counter()
    sequence = ''
    
    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    count_kmers_in_sequence(sequence, k, kmer_counts)
                    sequence = ''
            else:
                sequence += line
        if sequence:
            count_kmers_in_sequence(sequence, k, kmer_counts)
    return kmer_counts

def count_kmers_in_sequence(sequence, k, kmer_counts):    
    ascii_sequence = sequence.encode('ascii')
    kmer_encoded = 0
    mask = (1 << (2 * k)) - 1
    kmer_size = 0  

    for ascii_char in ascii_sequence:
        code = encode_array[ascii_char]
        if code == -1:
            kmer_size = max(kmer_size - 1, 0)
            continue

        kmer_encoded = ((kmer_encoded << 2) | code) & mask
        if kmer_size < k:
            kmer_size += 1
        if kmer_size == k:
            kmer_counts[kmer_encoded] += 1
    
def generate_kmer(kmer_encoded, k):
    kmer = [''] * k
    offset = k-1

    for i in range(k):
        kmer[offset-i] = decode_list[kmer_encoded & 0b11]
        kmer_encoded >>= 2
    return ''.join(kmer)

def get_top_kmers(kmer_counts, k):
    counts_to_kmers = defaultdict(list)
    for kmer, count in kmer_counts.items():
        counts_to_kmers[count].append(kmer)

    sorted_counts = sorted(counts_to_kmers.keys(), reverse=True)

    result = []
    total_kmers = 0

    for count in sorted_counts:
        kmers = counts_to_kmers[count]
        kmers.sort()
        result.extend( (generate_kmer(kmer, k), count) for kmer in kmers )
        total_kmers += len(kmers)
        if total_kmers >= 100:
            break

    return result

def write_output(kmer_counts, output_file, k):
    top_100 = get_top_kmers(kmer_counts, k)   
    with open(output_file, 'w', newline='') as file:
        # lineterminator를 빈 문자열로 설정
        writer = csv.writer(file, lineterminator='\n')
        writer.writerows(top_100)

if __name__ == '__main__':
    main()
