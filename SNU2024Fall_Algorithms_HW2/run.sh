#!/bin/bash

# 사용법 체크
if [ "$#" -ne 3 ]; then
    echo "Usage: sh run.sh [k] [genome1.fasta] [genome2.fasta]"
    exit 1
fi

# 입력 인자 받기
k=$1
genome1=$2
genome2=$3

# 유전체 처리 함수 정의
process_genome() {
    local genome=$1
    local kmer_output="${genome%.fasta}_filtered_kmers.txt"

    python cnt_kmer.py "$k" "$genome" | python filter_kmers.py "$k" > "$kmer_output"
}

# 백그라운드로 유전체 처리 시작
process_genome "$genome1" &
pid1=$!
process_genome "$genome2" &
pid2=$!

# 두 프로세스가 끝날 때까지 대기
wait $pid1
wait $pid2

# compute_lcs.py 실행 및 에러 핸들링
python compute_lcs.py "$genome1" "$genome2" "${genome1%.fasta}_filtered_kmers.txt" "${genome2%.fasta}_filtered_kmers.txt" "$k"