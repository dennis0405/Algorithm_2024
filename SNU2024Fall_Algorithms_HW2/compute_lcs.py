import sys

def read_filtered_kmers(file_path):
    with open(file_path, 'r') as f:
        return set(line.strip() for line in f if line.strip())

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        sequence = []
        for line in file:
            if line.startswith(">"):  # Skip headers
                continue
            sequence.append(line.strip().upper())
    return ''.join(sequence)

def find_kmer_positions(sequence, kmer_set, k):
    positions = []
    seq_length = len(sequence)
    for i in range(seq_length - k + 1):
        kmer = sequence[i:i+k] # 원본 sequence의 kmer
        
        if 'N' in kmer:
            if i+k < seq_length and kmer[0] != 'N':
                kmer_list = list(kmer)
                kmer_list.remove('N')
                kmer_list.append(sequence[i+k])
                kmer = ''.join(kmer_list)
        
        if kmer in kmer_set:
            position = i + 1
            positions.append((kmer, position))
    return positions

def compute_lcs(positions1, positions2):
    kmers1, positions_seq1 = zip(*positions1)
    kmers2, positions_seq2 = zip(*positions2)
    
    # k-mer 시퀀스를 기반으로 LCS 계산
    lcs, lcs_length = lcs_indices_kmers(kmers1, kmers2)

    if lcs_length == 0:
        return lcs_length, [], [], []

    # LCS에 해당하는 k-mer와 위치 가져오기
    lcs_kmers = []
    lcs_positions1 = []
    lcs_positions2 = []

    for i, j in lcs:
        lcs_kmers.append(kmers1[i])
        lcs_positions1.append(positions_seq1[i])
        lcs_positions2.append(positions_seq2[j])

    return lcs_length, lcs_kmers, lcs_positions1, lcs_positions2

def lcs_indices_kmers(X, Y):
    m, n = len(X), len(Y)
    prev_dp = [0]*(n+1)
    dp = [0]*(n+1)
    backtrack = [[0b00]*(n+1) for _ in range(m+1)]

    for i in range(1, m+1):
        for j in range(1, n+1):
            if X[i-1] == Y[j-1]:
                dp[j] = prev_dp[j-1]+1
                backtrack[i][j] = 0b11 # 대각선
            else:
                if prev_dp[j] >= dp[j-1]:
                    dp[j] = prev_dp[j]
                    backtrack[i][j] = 0b10 # 위로
                else:
                    dp[j] = dp[j-1] # 왼쪽은 0으로
        prev_dp = dp[:]

    lcs_length = dp[n]
    
    # LCS 추적
    lcs = []
    i, j = m, n
    while i > 0 and j > 0:
        if backtrack[i][j] == 0b11:
            lcs.append((i-1, j-1))  # 인덱스는 0-based
            i -= 1
            j -= 1
        elif backtrack[i][j] == 0b10:
            i -= 1
        else:
            j -= 1
    lcs.reverse()
    return lcs, lcs_length

def main():
    if len(sys.argv) != 6:
        print("Usage: python compute_lcs.py <genome1_fasta> <genome2_fasta> <filtered_kmer_file1> <filtered_kmer_file2> <k>")
        sys.exit(1)
    
    genome1_fasta = sys.argv[1]
    genome2_fasta = sys.argv[2]
    filter1 = sys.argv[3]
    filter2 = sys.argv[4]
    k = int(sys.argv[5])
    
    output_file_LCS = f"202116621_{k}_{genome1_fasta}_{genome2_fasta}_LCS.txt"
    output_file_genome1_LCS = f"202116621_{k}_{genome1_fasta}_LCS_positions.csv"
    output_file_genome2_LCS = f"202116621_{k}_{genome2_fasta}_LCS_positions.csv"

    # 필터링된 k-mer 불러오기
    kmer_set1 = read_filtered_kmers(filter1)
    kmer_set2 = read_filtered_kmers(filter2)
    kmer_set = kmer_set1.intersection(kmer_set2)

    # 원본 게놈 시퀀스 불러오기
    sequence1 = read_fasta(genome1_fasta)
    sequence2 = read_fasta(genome2_fasta)

    # k-mer 위치 찾기
    positions1 = find_kmer_positions(sequence1, kmer_set, k)
    positions2 = find_kmer_positions(sequence2, kmer_set, k)

    if not positions1 or not positions2:
        print("No k-mers left after filtering.")
        sys.exit(1)

    # LCS 계산
    lcs_length, lcs_kmers, lcs_positions1, lcs_positions2 = compute_lcs(positions1, positions2)

    if lcs_length == 0:
        print("No LCS found between the two genomes.")
        sys.exit(1)

    with open(output_file_LCS, 'w') as f:
        output = "-".join(lcs_kmers)
        f.write(output)

    with open(output_file_genome1_LCS, 'w') as f:
        f.write("LCS k-mer, Position\n")
        for kmer, position in zip(lcs_kmers, lcs_positions1):
            output = f"{kmer},{position}\n"
            f.write(output)
    
    with open(output_file_genome2_LCS, 'w') as f:
        f.write("LCS k-mer, Position\n")
        for kmer, position in zip(lcs_kmers, lcs_positions2):
            output = f"{kmer},{position}\n"
            f.write(output)
    
if __name__ == "__main__":
    main()
