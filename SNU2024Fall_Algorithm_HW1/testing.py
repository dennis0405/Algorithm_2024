import time

import os

# 스크립트의 디렉토리 경로 얻기
script_dir = os.path.dirname(os.path.abspath(__file__))

# 파일 경로 생성
file_path = os.path.join(script_dir, "Nonpathogenic_Escherichia coli ATCC 25922.fna")


# 방법 1: 일반적인 파일 읽기
start_time = time.time()
with open(file_path, 'r') as file:
    data = file.read()
print("일반 파일 읽기 시간:", time.time() - start_time)

# 방법 2: mmap 사용
import mmap

start_time = time.time()
with open(file_path, 'r') as file:
    with mmap.mmap(file.fileno(), length=0, access=mmap.ACCESS_READ) as mmapped_file:
        data = mmapped_file.read()
print("mmap 읽기 시간:", time.time() - start_time)
