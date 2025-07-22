import os

import pytest

from cosap._utils import get_bam_index_path, MultipleFilesFoundError


BAM_FILE = "alignment.bam"
BAI_FILE = "alignment.bai"
BAM_BAI_FILE = "alignment.bam.bai"


class TestGetBamIndexPath:
    def test_bam_invalid_extension(self, tmp_path):
        bai_path = os.path.join(tmp_path, BAI_FILE)

        with pytest.raises(ValueError, match="[Ee]xtension"):
            _ = get_bam_index_path(bai_path)


    def test_no_bam_file(self, tmp_path):
        bam_path = os.path.join(tmp_path, BAM_FILE)

        with pytest.raises(FileNotFoundError, match="BAM file"):
            _ = get_bam_index_path(bam_path)
    
    
    def test_no_bam_index_file(self, tmp_path):
        bam_path = os.path.join(tmp_path, BAM_FILE)
        open(bam_path, "w").write("alignment contents")

        with pytest.raises(FileNotFoundError, match="BAM index file"):
            _ = get_bam_index_path(bam_path)
    

    def test_bai_file_found(self, tmp_path):
        bam_path = os.path.join(tmp_path, BAM_FILE)
        bai_path = os.path.join(tmp_path, BAI_FILE)

        open(bam_path, "w").write("alignment contents")
        open(bai_path, "w").write("alignment index contents")

        assert get_bam_index_path(bam_path) == bai_path
    

    def test_bam_bai_file_found(self, tmp_path):
        bam_path = os.path.join(tmp_path, BAM_FILE)
        bam_bai_path = os.path.join(tmp_path, BAM_BAI_FILE)

        open(bam_path, "w").write("alignment contents")
        open(bam_bai_path, "w").write("alignment index contents")

        assert get_bam_index_path(bam_path) == bam_bai_path
    

    def test_multiple_index_files_found(self, tmp_path):
        bam_path = os.path.join(tmp_path, BAM_FILE)
        bai_path = os.path.join(tmp_path, BAI_FILE)
        bam_bai_path = os.path.join(tmp_path, BAM_BAI_FILE)

        open(bam_path, "w").write("alignment contents")
        open(bai_path, "w").write("alignment index contents")
        open(bam_bai_path, "w").write("alignment index contents")

        with pytest.raises(MultipleFilesFoundError):
            _ = get_bam_index_path(bam_path)
