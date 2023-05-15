from Bio import SeqIO

def extract_sequence(sequence, start, end):
    """
    从给定的序列中提取指定位置的子序列
    :param sequence: str, 输入的序列
    :param start: int, 起始位置
    :param end: int, 终止位置
    :return: str, 提取的子序列
    """
    return sequence[start-1:end]

# 读取fasta文件
fasta_file = "input.fasta"
fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

# 读取bed文件
bed_file = "input.bed"
with open(bed_file, "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        name = fields[3]
        
        # 提取子序列
        for fasta in fasta_sequences:
            if fasta.id == chrom:
                sequence = str(fasta.seq)
                subsequence = extract_sequence(sequence, start, end)
                
                # 将子序列写入输出文件
                output_file = name + ".fasta"
                with open(output_file, "w") as out:
                    out.write(">" + name + "\n")
                    out.write(subsequence + "\n")
                break

# 关闭fasta文件
fasta_sequences.close()
