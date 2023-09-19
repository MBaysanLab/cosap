import pytest
import os
from hashlib import sha1

from cosap._config import AppConfig

file_sizes = {
    "1000G_phase1.snps.high_confidence.hg38.vcf.gz": 1888262073,
    "1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi": 2128536,
    "Homo_sapiens_assembly38.1.bt2": 1019133809,
    "Homo_sapiens_assembly38.2.bt2": 760863372,
    "Homo_sapiens_assembly38.3.bt2": 40751,
    "Homo_sapiens_assembly38.4.bt2": 760863367,
    "Homo_sapiens_assembly38.dbsnp138.elsites": 1372365563,
    "Homo_sapiens_assembly38.dbsnp138.vcf": 10950827213,
    "Homo_sapiens_assembly38.dbsnp138.vcf.idx": 12480412,
    "Homo_sapiens_assembly38.dict": 581712,
    "Homo_sapiens_assembly38.elfasta": 3217493596,
    "Homo_sapiens_assembly38.fasta": 3249912778,
    "Homo_sapiens_assembly38.fasta.alt": 487553,
    "Homo_sapiens_assembly38.fasta.amb": 20199,
    "Homo_sapiens_assembly38.fasta.ann": 455474,
    "Homo_sapiens_assembly38.fasta.bwt": 3217347004,
    "Homo_sapiens_assembly38.fasta.bwt.2bit.64": 10456377594,
    "Homo_sapiens_assembly38.fasta.fai": 160928,
    "Homo_sapiens_assembly38.fasta.pac": 804336731,
    "Homo_sapiens_assembly38.fasta.sa": 1608673512,
    "Homo_sapiens_assembly38.rev.1.bt2": 1019133809,
    "Homo_sapiens_assembly38.rev.2.bt2": 760863372,
    "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz": 20685880,
    "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi": 1500013,
}
file_sha1_hashes = {
    "1000G_phase1.snps.high_confidence.hg38.vcf.gz": "5572a3946857ef76674ecde9a7959fb952167d44",
    "1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi": "f559a9d05a1990f2ef401738d35508a95fa278b3",
    "Homo_sapiens_assembly38.1.bt2": "d47aff77fdbc17505aa9aeb01ca44d71e7e57b26",
    "Homo_sapiens_assembly38.2.bt2": "fdb39086767f8a37807d80e706d66417a19367b1",
    "Homo_sapiens_assembly38.3.bt2": "3a740de75152e0cc8ebad656151a09660b26952a",
    "Homo_sapiens_assembly38.4.bt2": "d6c96e29ec2dfd2a1088cb24ee0ddbabcb8cd8fb",
    "Homo_sapiens_assembly38.dbsnp138.elsites": "2992cdb9ec2d34ae6087c06eb43a268d7518bde8",
    "Homo_sapiens_assembly38.dbsnp138.vcf": "6a6ba06158f69b82ce87c384bb638c782fd34c0b",
    "Homo_sapiens_assembly38.dbsnp138.vcf.idx": "9cadfd70f5fb1b215caa6074bcd52fe4c1aa9a24",
    "Homo_sapiens_assembly38.dict": "491a8fe834e818728105814df1870f0c8528e9bc",
    "Homo_sapiens_assembly38.elfasta": "081ba2ce32f3e9e49ad43d7ad79df0acd35a9668",
    "Homo_sapiens_assembly38.fasta": "2b9475654e19946c2067ee799f7a5b7115764b31",
    "Homo_sapiens_assembly38.fasta.alt": "7cbcd5f14b2f487f0167e0648152ed3a35dcacb6",
    "Homo_sapiens_assembly38.fasta.amb": "0854fc8106e40b5c8eb8288142444e59dc77eaa8",
    "Homo_sapiens_assembly38.fasta.ann": "3e9c41c6b5a7778ee6f3e8c17468fb6795bb2d93",
    "Homo_sapiens_assembly38.fasta.bwt": "d5eb32215cf3f00329f1c10325d1e4d05491ddc5",
    "Homo_sapiens_assembly38.fasta.bwt.2bit.64": "5603a2c3b39bb1f36102688a3db61d26c1bfc907",
    "Homo_sapiens_assembly38.fasta.fai": "521fa7ba606ad4cfbbf8fc6ac5ac95cccfcaf7a3",
    "Homo_sapiens_assembly38.fasta.pac": "18a8066606baa094479c94647d7556ae031edbce",
    "Homo_sapiens_assembly38.fasta.sa": "f90113ba3fa6ffe304e65f8d5b1276995bd4799c",
    "Homo_sapiens_assembly38.rev.1.bt2": "5d65af17c167d98e019ca4c2d8206860a5411c64",
    "Homo_sapiens_assembly38.rev.2.bt2": "c898c2158f6b6a4ec79923682e87d006d09a9e4c",
    "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz": "da89868391799231f9437b096a2afda6bd6e86ad",
    "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi": "a7d87f95887d4ceef98bb3b48a16e083665fdefe",
}

@pytest.mark.library
@pytest.mark.minute_long
def test_all_library_files():
    library_files = file_sizes.keys()
    
    for library_file in library_files:
        library_file_path = os.path.join(AppConfig.LIBRARY_PATH, library_file)
        assert os.path.isfile(library_file_path)
        assert os.path.getsize(library_file_path) == file_sizes[library_file]
        assert _hash_large_file_sha1(library_file_path) == file_sha1_hashes[library_file]


def _hash_large_file_sha1(file_path: str, block_size: int = 2**20) -> str:
    sha1_hasher = sha1()
    file = open(file_path, "rb")

    block = file.read(block_size)
    while block:
        sha1_hasher.update(block)
        block = file.read(block_size)
    
    return sha1_hasher.hexdigest()
