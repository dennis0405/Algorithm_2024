import sys
from collections import defaultdict

# k=10, Nonpathogenic 으로 했을 때 s (텍스트 파일 삭제)
# k=12, Pathogenic 으로 했을 때, 1m 10.345s 

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
    kmer_counts = defaultdict(lambda: ["", 0])
    current_sequence = ''
    
    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # 맨 첫 주석은 건너뜀, 현재 시퀀스가 있으면 처리
                if current_sequence:
                    count_kmers_in_sequence(current_sequence, k, kmer_counts, encode_array)
                    current_sequence = ''
            else:
                sequence_line = ''.join([c for c in line if c in 'ATCG'])
                current_sequence += sequence_line
        # 마지막 시퀀스 처리
        if current_sequence:
            count_kmers_in_sequence(current_sequence, k, kmer_counts, encode_array)
    return kmer_counts

def count_kmers_in_sequence(sequence, k, kmer_counts, encode_array):
    ascii_sequence = sequence.encode('ascii')
    kmer_encoded = 0  # encode된 DNA kmer
    mask = (1 << (2 * k)) - 1

    kmer = sequence[0:k]
    kmer_ascii = ascii_sequence[0:k]
    for ascii in kmer_ascii:
            kmer_encoded = (kmer_encoded << 2) | encode_array[ascii]
    kmer_counts[kmer_encoded][0] = kmer
    kmer_counts[kmer_encoded][1] += 1

    for i in range(1, len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmer_ascii = ascii_sequence[i+k-1]
        kmer_encoded = ((kmer_encoded << 2) | encode_array[kmer_ascii]) & mask
    
        kmer_counts[kmer_encoded][0] = kmer
        kmer_counts[kmer_encoded][1] += 1

    return kmer_counts

def generate_kmer(kmer_encoded, k):
    kmer = [''] * k
    offset = k-1

    for i in range(k):
        kmer[offset-i] = decode_list[kmer_encoded & 0b11]
        kmer_encoded >>= 2
    return ''.join(kmer)

def write_output(kmer_counts, output_file, k):
    max_kmer_count = (1 << (2 * k))
    # 입출력 개선
    buffer = []

    # kmer_int를 0부터 max_kmer_int-2까지 순회
    for i in range(max_kmer_count - 1):
        kmer, count = kmer_counts.get(i, ("", 0))

        if kmer == "":
            kmer = generate_kmer(i, k)
        
        buffer.append(f"{kmer},{count}\n")
    
    with open(output_file, 'w') as file:   
        file.writelines(buffer)

        kmer, count = kmer_counts.get(max_kmer_count-1, ("", 0))

        if kmer == "":
            kmer = generate_kmer(max_kmer_count-1, k)
        
        file.write(f"{kmer},{count}")
            
if __name__ == '__main__':
    main()
