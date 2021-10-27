import os


rule mark_dup:
    input:
        bam=lambda wildcards: os.path.normpath(
            os.path.join(
                config["workdir"],
                config["mapping"][wildcards.identification]["output"],
            )
        ),
    output:
        mdup_bam=os.path.normpath(
            os.path.join(config["workdir"], "{identification}_mdup.bam")
        ),
    log:
        os.path.normpath(
            os.path.join(config["workdir"], "{identification}_picard_markdup_log.txt")
        ),
    run:
        os.system(
            f"picard MarkDuplicates I={input.bam} O={output.mdup_bam} METRICS_FILE={wildcards.identification}_mdupmetrics.out"
        )


rule gatk_base_cal:
    input:
        bam=os.path.normpath(
            os.path.join(config["workdir"], "{identification}_mdup.bam")
        ),
        dbsnp=os.path.join(
            config["library_path"], "Homo_sapiens_assembly38.dbsnp138.vcf"
        ),
        mills_and_1000g=os.path.join(
            config["library_path"], "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        ),
        _1000g=os.path.join(
            config["library_path"], "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
        ),
        ref_genome=os.path.join(config["library_path"], "Homo_sapiens_assembly38.fasta"),
    output:
        calibration_table=os.path.normpath(
            os.path.join(config["workdir"], "{identification}_calibration_table.bam")
        ),
        calibrated_bam=os.path.normpath(
            os.path.join(config["workdir"], "{identification}_calibrated.bam")
        ),
    log:
        os.path.normpath(
            os.path.join(config["workdir"], "{identification}_gatk_base_cal_log.txt")
        ),
    run:
        os.system(
            f"""gatk BaseRecalibrator -R {input.ref_genome} -I {input.bam} \
        --known-sites {input.dbsnp} --known-sites {input.mills_and_1000g} --known-sites {input._1000g} -O {output.calibration_table}"""
        )
        os.system(
            f"gatk ApplyBQSR -R {input.ref_genome} -I {input.bam} --bqsr-recal-file {output.calibration_table} -O {output.calibrated_bam}"
        )
