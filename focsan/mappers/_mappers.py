import glob
import itertools as it
from abc import ABC, abstractmethod
from parser import Parser
from subprocess import run
from typing import Dict, List

from .._library_paths import LibraryPaths
from .._pipeline_config import PipelineConfig
from .._utils import join_paths


class _Mapper(ABC):
    @abstractmethod
    def map(self):
        pass


class _Mappable:
    RG_ID = "RG_ID"
    RG_SM = "RG_SM"
    RG_LB = "RG_LB"
    RG_PL = "RG_PL"
    RG_PU = "RG_PU"

    @classmethod
    def _list_fastq_files(cls, file_path: str) -> List:
        return glob.glob(os.path.join(file_path, "*fastq.gz"))

    @classmethod
    def _get_file_format_string(cls, sample_type: str, trimmed: bool) -> str:
        format_string = ""
        if trimmed:
            format_string = "{Trim}_"

        if sample_type == "tumor":
            format_string += (
                "{Sample_ID}_{Index}_{Lanes}_{Pair}_{Number_of_seq}.fastq.gz"
            )
        elif sample_type in ("germline", "normal"):
            format_string += (
                "{Sample_ID}_{Germline}_{Index}_{Lanes}_{Pair}_{Number_of_seq}.fastq.gz"
            )
        else:
            raise Exception(f"Unknown sample type: {sample_type}")

        return format_string

    @classmethod
    def _read_flowcell_info(self, file_path: str) -> str:
        with gzip.open(file_path) as fp:
            header = str(fp.readline())
        if header.startswith("@"):
            raise Exception(f"Missing header on file {file_path}.")

        return header.split(":")[2]

    @classmethod
    def _get_file_information(
        cls,
        pipeline_config: PipelineConfig,
    ) -> List:
        fastq_list = [
            join_paths(pipeline_config.FASTQ_DIR, filename)
            for filename in glob.glob(
                join_paths(pipeline_config.FASTQ_DIR, "*fastq.gz")
            )
        ]
        format_string = cls._get_file_format_string(
            pipeline_config.SAMPLE_TYPE, pipeline_config.FASTQ_TRIMMED
        )
        fastq_filename_parser = Parser(format_string)
        info_dict = defaultdict(dict)
        flowcell_info = cls._read_flowcell_info(fastq_list[0])
        for file_path in fastq_list:
            filename = os.path.basename(file_path)
            file_info = {
                **fastq_filename_parser.parse(filename).named,
                file_info["Pair"]: file_path,
                "Flowcell": flowcell_info,
            }
            info_dict[(file_info["Lanes"], file_info["Number_of_seq"])].update(
                file_info
            )

        return list(info_dict.values())

    @classmethod
    def _create_output_filename(
        cls, file_info: Dict, pipeline_config: PipelineConfig
    ) -> str:
        return f"""{pipeline_config.MAPPER_TYPE}_{fastq_info["Sample_ID"]}_{fastq_info["Index"]}_{fastq_info["Lanes"]}_{fastq_info["Number_of_seq"]}.bam"""

    @classmethod
    def _get_rg_flags(cls, fastq_info: Dict, pipeline_config: PipelineConfig) -> Dict:
        flags = {
            cls.RG_ID: f"""{fastq_info["Flowcell"]}.{fastq_info["Lanes"][-1]}""",
            cls.RG_SM: fastq_info["Sample_ID"],
            cls.RG_LB: pipeline_config.PATIENT_ID,
            cls.RG_PL: "Illumina",
            cls.RG_PU: f"""{fastq_info["Flowcell"]}.{fastq_info["Index"]}.{fastq_info["Lanes"][-1]}""",
        }
        return flags


class BWAMapper(_Mapper, _Mappable):
    @classmethod
    def _create_read_group(
        cls, fastq_info: Dict, pipeline_config: PipelineConfig
    ) -> str:
        flags = cls._get_rg_flags(
            fastq_info=fastq_info, pipeline_config=pipeline_config
        )
        read_arguments = "".join(
            (
                '"',
                r"@RG\tID",
                flags[cls.RG_ID],
                r"\tSM:",
                flags[cls.RG_SM],
                r"\tLB",
                flags[cls.RG_LB],
                r"\tPL",
                flags[cls.RG_PL],
                r"\tPU",
                flags[cls.RG_PU],
                '"',
            )
        )
        return read_arguments

    @classmethod
    def _create_command(
        cls,
        pipeline_config: PipelineConfig,
        fastq_info: Dict,
        library_paths: LibraryPaths,
    ) -> List:
        output_filename = cls._create_output_filename(
            fastq_info=fastq_info, pipeline_config=pipeline_config
        )
        read_group = cls._create_read_group(
            fastq_info=fastq_info, pipeline_config=pipeline_config
        )
        command = [
            "bwa",
            "mem",
            "-t",
            pipeline_config.MAPPER_THREADS,
            "-R",
            read_group,
            library_paths.REF_DIR,
            # TODO: This needs to go somewhere else
            os.path.normpath("Bwa/Homo_sapiens_assembly38.fasta"),
            fastq_info["R1"],
            fastq_info["R2"],
            "|",
            "samtools",
            "view",
            "-@",
            pipeline_config.MAPPER_THREADS,
            "-bS",
            "-",
            ">",
            output_filename,
        ]
        return command

    @classmethod
    def map(cls, pipeline_config: PipelineConfig):
        library_paths = LibraryPaths()
        fastq_info_list = cls._get_file_information(pipeline_config=pipeline_config)

        for fastq_info in fastq_info_list:
            command = cls._create_command(
                pipeline_config=pipeline_config,
                fastq_info=fastq_info,
                library_paths=library_paths,
            )
            run(command, cwd=pipeline_config.FASTQ_DIR)


class Bowtie2Mapper(_Mapper, _Mappable):
    @classmethod
    def _create_read_group(
        cls, fastq_info: Dict, pipeline_config: PipelineConfig
    ) -> List:
        flags = cls._get_rg_flags(
            fastq_info=fastq_info, pipeline_config=pipeline_config
        )
        read_arguments = [
            "--rg-id",
            flags[cls.RG_ID],
            "--rg",
            f"SM:{flags[cls.RG_SM]}",
            "--rg",
            f"LB:{flags[cls.RG_LB]}",
            "--rg",
            f"PL:{flags[cls.RG_PL]}",
            "--rg",
            f"PU:{flags[cls.RG_PU]}",
        ]
        return read_arguments

    @classmethod
    def _create_command(
        cls,
        pipeline_config: PipelineConfig,
        fastq_info: Dict,
        library_paths: LibraryPaths,
    ) -> List:
        output_filename = cls._create_output_filename(
            fastq_info=fastq_info, pipeline_config=pipeline_config
        )
        read_group = cls._create_read_group(
            fastq_info=fastq_info, pipeline_config=pipeline_config
        )
        command = [
            "bowtie2",
            "-p",
            pipeline_config.MAPPER_THREADS,
            *read_group,
            "-x",
            library_paths.REF_DIR,
            # TODO: This needs to go somewhere else
            os.path.normpath("Bowtie2/Homo_sapiens_assembly38"),
            "-1",
            fastq_info["R1"],
            "-2",
            fastq_info["R2"],
            "|",
            "samtools",
            "view",
            "-@",
            pipeline_config.MAPPER_THREADS,
            "-bS",
            "-",
            ">",
            output_filename,
        ]
        return command

    @classmethod
    def map(cls, pipeline_config: PipelineConfig):
        library_paths = LibraryPaths()
        fastq_info_list = cls._get_file_information(pipeline_config=pipeline_config)

        for fastq_info in fastq_info_list:
            command = cls._create_command(
                pipeline_config=pipeline_config,
                fastq_info=fastq_info,
                library_paths=library_paths,
            )
            run(command, cwd=pipeline_config.FASTQ_DIR)


class NovoalignMapper(_Mapper, _Mappable):
    @classmethod
    def _create_read_group(
        cls, fastq_info: Dict, pipeline_config: PipelineConfig
    ) -> str:
        flags = cls._get_rg_flags(
            fastq_info=fastq_info, pipeline_config=pipeline_config
        )
        read_arguments = "".join(
            (
                '"',
                r"@RG\tID",
                flags[cls.RG_ID],
                r"\tSM:",
                flags[cls.RG_SM],
                r"\tLB",
                flags[cls.RG_LB],
                r"\tPL",
                flags[cls.RG_PL],
                r"\tPU",
                flags[cls.RG_PU],
                '"',
            )
        )
        return read_arguments

    @classmethod
    def _create_command(
        cls,
        pipeline_config: PipelineConfig,
        fastq_info: Dict,
        library_paths: LibraryPaths,
    ) -> List:
        output_filename = cls._create_output_filename(
            fastq_info=fastq_info, pipeline_config=pipeline_config
        )
        read_group = cls._create_read_group(
            fastq_info=fastq_info, pipeline_config=pipeline_config
        )
        stats_filename = f"{os.path.splitext(output_filename)[0]}_stats.txt"
        command = [
            "novoalign",
            "-k",
            "-d",
            library_paths.REF_DIR,
            os.path.normpath("Bowtie2/Homo_sapiens_assembly38"),
            "-f",
            fastq_info["R1"],
            fastq_info["R2"],
            "-a",
            "-c",
            pipeline_config.MAPPER_THREADS,
            "-o",
            "SAM",
            *read_group,
            "2>",
            stats_filename,
            "|",
            "samtools",
            "view",
            "-@",
            pipeline_config.MAPPER_THREADS,
            "-bS",
            "-",
            ">",
            output_filename,
        ]
        return command

    @classmethod
    def map(cls, pipeline_config: PipelineConfig):
        library_paths = LibraryPaths()
        fastq_info_list = cls._get_file_information(pipeline_config=pipeline_config)

        for fastq_info in fastq_info_list:
            command = cls._create_command(
                pipeline_config=pipeline_config,
                fastq_info=fastq_info,
                library_paths=library_paths,
            )
            run(command, cwd=pipeline_config.FASTQ_DIR)
