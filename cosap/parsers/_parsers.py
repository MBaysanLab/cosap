from abc import ABC, abstractmethod


class _Parser(ABC):
    @abstractmethod
    def parse(self):
        pass


class _Parsable:
    pass


class FastqParser(_Parser, _Parsable):
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
            file_info = fastq_filename_parser.parse(filename).named
            file_info.update(
                {
                    file_info["Pair"]: file_path,
                    "Flowcell": flowcell_info,
                }
            )
            info_dict[(file_info["Lanes"], file_info["Number_of_seq"])].update(
                file_info
            )

        return list(info_dict.values())

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

    @classmethod
    def _create_output_filename(
        cls, file_info: Dict, pipeline_config: PipelineConfig
    ) -> str:
        return f"""{pipeline_config.MAPPER_TYPE}_{file_info["Sample_ID"]}_{file_info["Index"]}_{file_info["Lanes"]}_{file_info["Number_of_seq"]}.bam"""

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
