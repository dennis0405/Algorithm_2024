import sys
from collections import defaultdict


# k=10, Nonpathogenic 으로 했을 때 6.263s (텍스트 파일 삭제)
# k=12, Pathogenic 으로 했을 때, 51.362s 

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
    
    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    count_kmers_in_sequence(current_sequence, k, kmer_counts)
                    current_sequence = ''
            else:
                sequence_line = ''.join([c for c in line if c in 'ATCG'])
                current_sequence += sequence_line
        if current_sequence:
            count_kmers_in_sequence(current_sequence, k, kmer_counts)
    return kmer_counts

def count_kmers_in_sequence(sequence, k, kmer_counts):
    
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmer_counts[kmer] += 1

    return kmer_counts

def encode_kmer(kmer):
    # A: 00, C: 01, G: 10, T: 11 의 2비트 인코딩을 위한 딕셔너리
    encode_dict = {
        'A': 0b00,
        'C': 0b01,
        'G': 0b10,
        'T': 0b11
    }
    
    encoded_value = 0
    for char in kmer:
        encoded_value = (encoded_value << 2) | encode_dict[char]  # 왼쪽으로 2비트 시프트 후 인코딩 값 추가
    
    return encoded_value

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

    # Counting Sort를 위한 배열 초기화
    count_array = [("", 0)] * max_kmer_count

    # kmer_counts의 값을 count_array에 복사
    for kmer, count in kmer_counts.items():
        encoded_kmer = encode_kmer(kmer)
        count_array[encoded_kmer] = (kmer, count)
    
    # kmer_int를 0부터 max_kmer_int-1까지 순회
    for i in range(max_kmer_count - 1):
        kmer_count = count_array[i]
        if kmer_count[0] == '':
            kmer = generate_kmer(i, k)
            buffer.append(f"{kmer},{kmer_count[1]}\n")

        buffer.append(f"{kmer_count[0]},{kmer_count[1]}\n")
    
    with open(output_file, 'w') as file:   
        file.writelines(buffer)
        kmer_count = count_array[max_kmer_count - 1]
        if kmer_count[0] == '':
            kmer = generate_kmer(i, k)
            file.write(f"{kmer},{kmer_count[1]}\n")

        file.write(f"{kmer_count[0]},{kmer_count[1]}")
            
if __name__ == '__main__':
    main()
