import sys
from collections import defaultdict
from itertools import product

# k=10, Nonpathogenic 으로 했을 때 s (텍스트 파일 삭제)
# k=12, Pathogenic 으로 했을 때, 17.133s 

encode_array = [-1] * 128     
encode_array[65] = 0b00  
encode_array[67] = 0b01  
encode_array[71] = 0b10  
encode_array[84] = 0b11  

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
    kmer_counts = defaultdict(int)
    current_sequence = ''
    
    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    count_kmers_in_sequence(current_sequence, k, kmer_counts, encode_array)
                    current_sequence = ''
            else:
                current_sequence += line
        if current_sequence:
            count_kmers_in_sequence(current_sequence, k, kmer_counts, encode_array)
    return kmer_counts

def count_kmers_in_sequence(sequence, k, kmer_counts, encode_array):
    s_ascii = sequence.encode('ascii')
    
    index = 0
    window = 0
    mask = (1 << (2 * k)) - 1  
    kmer_encoded = 0 

    while True:
        code = encode_array[s_ascii[index]]
        if(code == -1):
            index += 1
            continue
        kmer_encoded = (kmer_encoded << 2) | code
        window += 1
        index += 1
        if (window == k):    
            kmer_counts[kmer_encoded] += 1
            break
        
    for i in range(index, len(sequence)):
        code = encode_array[s_ascii[i]]
        if(code == -1):
            continue
        kmer_encoded = ((kmer_encoded << 2) | code) & mask
        kmer_counts[kmer_encoded] += 1

def generate_all_kmers(k):
    DNA = ['A', 'C', 'G', 'T']
    return [''.join(p) for p in product(DNA, repeat=k)]

def generate_kmer(kmer_encoded, k):
    kmer = [''] * k
    offset = k-1

    for i in range(k):
        kmer[offset-i] = decode_list[kmer_encoded & 0b11]
        kmer_encoded >>= 2
    return ''.join(kmer)

def write_output(kmer_counts, output_file, k):
    max_kmer_count = 1 << (2 * k)
    kmers = generate_all_kmers(k)
    
    count_array = [0] * max_kmer_count
    for kmer_encoded, count in kmer_counts.items():
        count_array[kmer_encoded] = count
    
    with open(output_file, 'w') as file:
        for kmer_encoded, kmer in enumerate(kmers):
            count = count_array[kmer_encoded]
            if kmer_encoded < max_kmer_count - 1:
                file.write(f"{kmer},{count}\n")
            else:
                file.write(f"{kmer},{count}")
            
if __name__ == '__main__':
    main()
