rule varscan:
    input:
        germline_bam=os.path.normpath(
            os.path.join(config["workdir"], "{germline_id}_calibrated.bam")
        ),
        tumor_bam=os.path.normpath(
            os.path.join(config["workdir"], "{tumor_id}_calibrated.bam")
        ),
        ref_genome=os.path.join(config["library_path"], "Homo_sapiens_assembly38.fasta"),
    output:
        snp_vcf=os.path.normpath(
            os.path.join(config["workdir"], "snp_{tumor_id}_{germline_id}_varscan.vcf")
        ),
    params:
        germline_sample_name=lambda wildcards: config["mapping"][
            wildcards.germline_id
        ]["params"]["read_groups"]["ID"],
    log:
        os.path.normpath(
            os.path.join(config["workdir"], "{tumor_id}_{germline_id}_varscan_log.txt")
        ),
    run:
        command_mpile_ = f"""samtools mpileup -f {input.ref_genome} {input.germline_bam} {input.tumor_bam} \
                        | varscan somatic --output-snp {output.snp_vcf} --mpileup 1 --output-vcf nthreads=8"""

        command_process_somatic = f"""varscan processSomatic {output.snp_vcf}"""

        os.system(command_mpile_)
        os.system(command_process_somatic)
