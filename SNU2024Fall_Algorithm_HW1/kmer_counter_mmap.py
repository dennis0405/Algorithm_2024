import sys
from collections import defaultdict
import mmap

encode_array = [-1] * 128      # 모든 ASCII 문자를 포함, 'ACGT'중 하나가 아닐 경우, -1로 init
encode_array[65] = 0b00  # encode_array[ord('A')] = 0b00
encode_array[67] = 0b01  # encode_array[ord('C')] = 0b01
encode_array[71] = 0b10  # encode_array[ord('G')] = 0b10
encode_array[84] = 0b11  # encode_array[ord('T')] = 0b11

decode_list = ['A', 'C', 'G', 'T']

def main():
    if len(sys.argv) != 3:
        print("Usage: python kmer_counter.py [k-mer] [input_file]")
        sys.exit(1)

    k = int(sys.argv[1])
    input_file = sys.argv[2]
    
    student_id = '202116621'
    output_file = f"{student_id}.txt"
    
    # 주어진 k값과 input_file을 바탕으로 kmer을 count (dict 형태)
    kmer_counts = process_file(k, input_file)
    
    # 결과를 출력 파일에 작성
    write_output(kmer_counts, output_file, k)


def process_file(k, input_file):
    kmer_counts = defaultdict(int)
    current_sequence = ''
    
    with open(input_file, 'r') as f:
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            # 굳이 문자열로 decode하지 않고, 바로 ascii 형태로 사용
            current_sequence = b''
            for line in iter(mm.readline, b""):
                line = line.strip()
                if line.startswith(b'>'):
                    if current_sequence:
                        count_kmers_in_sequence(current_sequence, k, kmer_counts, encode_array)
                        current_sequence = b''
                else:
                    current_sequence += line
            if current_sequence:
                count_kmers_in_sequence(current_sequence, k, kmer_counts, encode_array)
    return kmer_counts

def count_kmers_in_sequence(s, k, kmer_counts, encode_array):
    # encode_array에 바로 접근하기 위해, ascii code 형태로 encode된 byte array 사용
    index = 0
    window = 0
    mask = (1 << (2 * k)) - 1  # 0111...11 (1이 2k개 만큼 존재, 하위 2k bit만 사용하기 위함)
    kmer_encoded = 0  # encode된 DNA kmer

    while True:
        code = encode_array[s[index]]

        if(code == -1):
            index += 1
            continue
        
        kmer_encoded = (kmer_encoded << 2) | code
        window += 1
        index += 1

        if (window == k):    
            kmer_counts[kmer_encoded] += 1
            break
    
    for i in range(index, len(s)):
        code = encode_array[s[i]]

        if(code == -1):
            continue
        
        kmer_encoded = ((kmer_encoded << 2) | code) & mask
        kmer_counts[kmer_encoded] += 1


def decode_kmer(kmer_int, k):
    kmer = [''] * k
    for i in range(k - 1, -1, -1):
        bits = (kmer_int >> (2 * i)) & 0b11
        kmer[k-1-i] = decode_list[bits]
    return ''.join(kmer)


def write_output(kmer_counts, output_file, k):
    sorted_kmers = sorted(kmer_counts.items())
    with open(output_file, 'w') as f:
        for kmer_encoded, count in sorted_kmers:
            kmer = decode_kmer(kmer_encoded, k)
            f.write(f"{kmer},{count}\n")


if __name__ == '__main__':
    main()
