from typing import List, Dict
import os


def get_fastq_files(wildcards) -> List:
    fastq_files = list()
    for fastq in config["mapping"][wildcards.identification]["input"]:
        fastq_files.append(config["mapping"][wildcards.identification]["input"][fastq])

    return fastq_files


rule bwa:
    input:
        fastq=get_fastq_files,
        ref_genome=os.path.join(config["library_path"], "Homo_sapiens_assembly38.fasta"),
    params:
        read_groups=lambda wildcards: config["mapping"][wildcards.identification][
            "params"
        ]["read_groups"],
    output:
        bam=os.path.normpath(
            os.path.join(config["workdir"], "{identification}_bwa.bam")
        ),
    log:
        os.path.normpath(os.path.join(config["workdir"], "{identification}_log.txt")),
    run:
        read_group_str = (
            ' -R "@RG\\tID:'
            + params.read_groups["ID"]
            + "\\tSM:"
            + params.read_groups["SM"]
            + "\\tLB:"
            + params.read_groups["LB"]
            + "\\tPL:"
            + params.read_groups["PL"]
            + "\\tPU:"
            + params.read_groups["PU"]
            + '"'
        )

        command = f"bwa mem -t 8 {read_group_str} {input.ref_genome} {' '.join(input.fastq)} | samtools sort -@ 8 -o {output.bam} -"
        os.system(command)

