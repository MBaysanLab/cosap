"""Regression tests for the standalone VEP VCF to Arrow script."""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pyarrow.ipc as ipc


def _load_script_module(project_root):
    script_path = Path(project_root) / "scripts" / "vep_vcf_to_arrow.py"
    spec = spec_from_file_location("vep_vcf_to_arrow", script_path)
    module = module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_convert_excludes_csq_from_generic_info_columns(tmp_path, project_root):
    module = _load_script_module(project_root)

    input_path = tmp_path / "input.vcf"
    output_path = tmp_path / "output.arrow"
    input_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                (
                    "##INFO=<ID=CSQ,Number=.,Type=String,Description=\"Consequence"
                    " annotations from Ensembl VEP. Format: "
                    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|"
                    "BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|"
                    "Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|"
                    "STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|PICK\">"
                ),
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                (
                    "1\t12345\t.\tA\tG\t.\tPASS\t"
                    "DP=42;CSQ=G|missense_variant|MODERATE|GENE1|ENSG000001|"
                    "Transcript|ENST000001|protein_coding|2/5||ENST000001.2:c.123A>G|"
                    "ENSP000001.2:p.Lys41Arg|123|123|41|K/R|aAg/aGg|rs1||1||HGNC|1234|"
                    "YES|MANE1|1,G|synonymous_variant|LOW|GENE1|ENSG000001|Transcript|"
                    "ENST000002|protein_coding|2/5||ENST000002.2:c.123A>G|"
                    "ENSP000002.2:p.Lys41Lys|123|123|41|K/K|aAg/aAa|rs1||1||HGNC|1234|"
                    "YES|MANE1|0"
                ),
                "",
            ]
        ),
        encoding="utf-8",
    )

    module.convert(str(input_path), str(output_path), batch_size=1)

    with ipc.open_stream(output_path) as reader:
        table = reader.read_all()

    assert "info_dp" in table.column_names
    assert "info_csq" not in table.column_names
    assert table.num_rows == 2
    assert table.column("info_dp").to_pylist() == ["42", "42"]
