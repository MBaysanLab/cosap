from typing import List, Dict
import os

scattered_interval_files = [
    interval for interval in os.listdir(os.path.join(config["library_path"], "scattered_intervals"))
]

rule mutect:
    input:
        germline_bam=os.path.normpath(
            os.path.join(config["workdir"], "{germline_id}_calibrated.bam")
        ),
        tumor_bam=os.path.normpath(
            os.path.join(config["workdir"], "{tumor_id}_calibrated.bam")
        ),
        ref_genome=os.path.join(config["library_path"], "Homo_sapiens_assembly38.fasta"),
    output:
        vcf=os.path.normpath(
            os.path.join(config["workdir"], "mutect/scattered/{tumor_id}_{germline_id}_{interval}_mutect.vcf")
        )
    params:
        germline_sample_name=lambda wildcards: config["mapping"][
            wildcards.germline_id
        ]["params"]["read_groups"]["ID"],
    log:
        os.path.normpath(
            os.path.join(config["workdir"], "mutect/{tumor_id}_{germline_id}_{interval}_mutect_log.txt")
        ),
    run:
        interval_file = [interval for interval in scattered_interval_files if wildcards.interval in interval][0]
        interval_file = os.path.normpath(os.path.join(config["library_path"], "scattered_intervals", interval_file))
        
        command = f"""gatk Mutect2 -R {input.ref_genome} -I {input.tumor_bam}\
                             -I {input.germline_bam} -O {output.vcf} -normal {params.germline_sample_name}\
                             -intervals {interval_file} > {log}"""

        os.system(command)



def get_scattered_inputs(wildcards):
    vcf_list = []
    for interval in scattered_interval_files:
        interval = interval.split(".")[0].split("-")[0]
        vcf = os.path.normpath(os.path.join(config["workdir"], f"mutect/scattered/{wildcards.tumor_id}_{wildcards.germline_id}_{interval}_mutect.vcf"))
        vcf_list.append(vcf)
    return vcf_list


rule gather_scattered_calls:
    input:
        scattered_vcfs = get_scattered_inputs
    output:
        vcf=os.path.normpath(
            os.path.join(config["workdir"], "mutect/{tumor_id}_{germline_id}_merged.vcf")
        )
    log:
        os.path.normpath(
            os.path.join(config["workdir"], "mutect/{tumor_id}_{germline_id}_gather_vcfs_log.txt")
        ),
    run:
        inputs = [f"-I={vcf}" for vcf in input.scattered_vcfs]
        command = f"""picard GatherVcfs {inputs} -O {output.vcf}"""

rule select_variants:
    input:
        ref_genome=os.path.join(config["library_path"], "Homo_sapiens_assembly38.fasta"),
        vcf=rules.gather_scattered_calls.output.vcf
    output:
        snp_vcf=os.path.normpath(
            os.path.join(config["workdir"], "mutect/snp_{tumor_id}_{germline_id}_mutect.vcf")
        )
    log:
        os.path.normpath(
            os.path.join(config["workdir"], "mutect/{tumor_id}_{germline_id}_mutect_select_variants_log.txt")
        ),
    run:
        select_variants_command = f"""gatk SelectVariants -R {input.ref_genome} -V {input.vcf} --select-type-to-include SNP -O {output.snp_vcf} > {log}"""

        os.system(command)