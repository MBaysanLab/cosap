from typing import List, Dict
import os


def get_fastq_files(wildcards) -> List:
    fastq_files = list()
    for fastq in config["mapping"][wildcards.identification]["input"]:
        fastq_files.append(config["mapping"][wildcards.identification]["input"][fastq])

    return fastq_files


rule bowtie:
    input:
        fastq=get_fastq_files,
        ref_genome=os.path.join(config["library_path"], "Homo_sapiens_assembly38.fasta"),
    params:
        read_groups=lambda wildcards: config["mapping"][wildcards.identification][
            "params"
        ]["read_groups"],
    output:
        bam=os.path.normpath(
            os.path.join(config["workdir"], "{identification}_bowtie.bam")
        ),
    log:
        os.path.normpath(os.path.join(config["workdir"], "{identification}_log.txt")),
    run:
        read_group_str = (
            "--rg-id "
            + params.read_groups["ID"]
            + " --rg SM:"
            + params.read_groups["SM"]
            + " --rg LB:"
            + params.read_groups["LB"]
            + " --rg PL:"
            + params.read_groups["PL"]
            + " --rg PU:"
            + params.read_groups["PU"]
        )

        fastq_pairs = f"-1 {input.fastq[0]} -2 {input.fastq[1]}"

        bowtie_index = input.ref_genome.split(".")[0]

        command = f"bowtie2 --threads 8 -x {bowtie_index} {fastq_pairs} {read_group_str} | samtools sort -@ 8 -o {output.bam} -"
        os.system(command)
