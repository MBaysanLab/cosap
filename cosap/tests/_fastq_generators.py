import os

from cosap._library_paths import LibraryPaths


def create_test_fastqs_from_chr1_ref_with_1_snp(
    workdir,
):  # Returns SNP (CHROM, POS, REF, ALT)
    sample_name = "TEST_SAMPLE_"

    position_of_last_read_base = 0

    line_length = fasta_line_length(LibraryPaths().REF_FASTA)
    with open(LibraryPaths().REF_FASTA) as ref:
        ref.readline()  # Skip header line

        start_of_base_read_lines = ref.tell()
        while "N" in ref.readline():
            start_of_base_read_lines = ref.tell()
            position_of_last_read_base += line_length

        ref.seek(start_of_base_read_lines)
        seq = ""
        for _ in range(10):
            seq += ref.readline()
            position_of_last_read_base += line_length

    seq = seq.replace("\n", "")

    quality_chr = phred_quality_score_encoding(63)

    with open(os.path.join(workdir, f"{sample_name}N_1.fastq"), "w") as N_1:
        N_1.write(f"@{sample_name}N.1 1/1\n")
        N_1.write(f"{seq}\n")
        N_1.write(f"+\n")
        N_1.write(f"{quality_chr * len(seq)}\n")

    with open(os.path.join(workdir, f"{sample_name}N_2.fastq"), "w") as N_2:
        N_2.write(f"@{sample_name}N.1 1/2\n")
        N_2.write(f"{complement(reverse(seq))}\n")
        N_2.write(f"+\n")
        N_2.write(f"{quality_chr * len(seq)}\n")

    bad_seq = seq[:-1] + complement(seq[-1])

    with open(os.path.join(workdir, f"{sample_name}T_1.fastq"), "w") as T_1:
        T_1.write(f"@{sample_name}T.1 1/1\n")
        T_1.write(f"{bad_seq}\n")
        T_1.write(f"+\n")
        T_1.write(f"{quality_chr * len(bad_seq)}\n")

    with open(os.path.join(workdir, f"{sample_name}T_2.fastq"), "w") as T_2:
        T_2.write(f"@{sample_name}T.1 1/2\n")
        T_2.write(f"{complement(reverse(bad_seq))}\n")
        T_2.write(f"+\n")
        T_2.write(f"{quality_chr * len(bad_seq)}\n")

    return (
        (N_1.name, N_2.name),
        (T_1.name, T_2.name),
        ("chr1", position_of_last_read_base, seq[-1], bad_seq[-1]),
    )


def complement(seq: str) -> str:
    return (
        seq.replace("A", "1")
        .replace("T", "2")
        .replace("G", "3")
        .replace("C", "4")
        .replace("1", "T")
        .replace("2", "A")
        .replace("3", "C")
        .replace("4", "G")
    )


def reverse(seq: str) -> str:
    return seq[::-1]


def phred_quality_score_encoding(score: int) -> chr:
    return chr(score + 33)


def fasta_line_length(fasta_path) -> int:
    with open(fasta_path) as fasta:
        fasta.readline()
        return len(fasta.readline().replace("\n", ""))
