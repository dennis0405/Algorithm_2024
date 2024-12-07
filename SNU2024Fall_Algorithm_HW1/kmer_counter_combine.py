import sys
from collections import Counter

encode_array = [-1] * 128
encode_array[65] = 0b00  
encode_array[67] = 0b01  
encode_array[71] = 0b10  
encode_array[84] = 0b11 

decode_list = ['A', 'C', 'G', 'T']

def main():
    if len(sys.argv) != 3:
        print("Usage: python kmer_counter.py [k-mer] [input_file]")
        sys.exit(1)
    k = int(sys.argv[1])
    input_file = sys.argv[2]
    student_id = '202116621'
    output_file = f"{student_id}.txt"
    kmer_counts = process_file(k, input_file)
    write_output(kmer_counts, output_file, k)

def process_file(k, input_file):
    kmer_counts = Counter()
    sequence = []
    
    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    count_kmers_in_sequence(''.join(sequence), k, kmer_counts)
                    sequence = []
            else:
                sequence.append(line)
        if sequence:
            count_kmers_in_sequence(''.join(sequence), k, kmer_counts)
    return kmer_counts

def count_kmers_in_sequence(sequence, k, kmer_counts):    
    ascii_sequence = sequence.encode('ascii')
    window = 0
    kmer_encoded = 0
    mask = (1 << (2 * k)) - 1
    
    index = 0
    while True:
        code = encode_array[ascii_sequence[index]]
        if code == -1:
            index += 1
            continue

        kmer_encoded = (kmer_encoded << 2) | code
        window += 1
        index += 1

        if window == k:
            kmer_counts[kmer_encoded] += 1
            break  

    for i in range(index, len(sequence)):
        code = encode_array[ascii_sequence[i]]
        if code == -1:
            continue
    
        kmer_encoded = ((kmer_encoded << 2) | code) & mask
        kmer_counts[kmer_encoded] += 1
    
def generate_kmer(kmer_encoded, k):
    kmer = [''] * k
    offset = k-1

    for i in range(k):
        kmer[offset-i] = decode_list[kmer_encoded & 0b11]
        kmer_encoded >>= 2
    return ''.join(kmer)
    counts_to_kmers = defaultdict(list)
    for kmer, count in kmer_counts.items():
        counts_to_kmers[count].append(kmer)

    sorted_counts = sorted(counts_to_kmers.keys(), reverse=True)

    result = []
    total_kmers = 0

    for count in sorted_counts:
        kmers = counts_to_kmers[count]
        kmers.sort()
        result.extend( (kmer, count) for kmer in kmers )
        total_kmers += len(kmers)
        if total_kmers >= N:
            break

    return result

def write_output(kmer_counts, output_file, k):
    combined_list = []
    for kmer, count in kmer_counts.items():
        combined = (count << 32) - kmer 
        combined_list.append(combined)

    combined_list.sort(reverse=True)

    if len(combined_list) > 100:
        edge_combined = combined_list[99]
        edge_count = edge_combined >> 32
    else:
        edge_count = 0

    result = []
    for combined in combined_list:
        count = combined >> 32
        if count < edge_count:
            break
        kmer = - (combined - (count << 32))
        result.append((count, kmer))

    with open(output_file, 'w') as file:
        length = len(result) - 1
        i = 0
        for count, kmer_encoded in result:
            kmer = generate_kmer(kmer_encoded, k)
            line = f"{kmer},{count}"
            if i < length:
                file.write(line + "\n")
            else:
                file.write(line)
            i += 1
    
if __name__ == '__main__':
    main()
