import os
from tempfile import TemporaryDirectory

import pytest

from cosap import MemoryHandler, AppConfig
from cosap._utils import join_paths


# Prevent unexpected behaviour due to AppConfig being global
@pytest.fixture(autouse=True)
def manage_app_config():
    IN_MEMORY_MODE_backup = AppConfig.IN_MEMORY_MODE

    yield

    AppConfig.IN_MEMORY_MODE = IN_MEMORY_MODE_backup


@pytest.fixture
def memory_handler():
    AppConfig.IN_MEMORY_MODE = False
    return MemoryHandler()


@pytest.fixture
def memory_handler_in_memory():
    AppConfig.IN_MEMORY_MODE = True
    return MemoryHandler()


def test_init(memory_handler):
    assert os.path.isdir(memory_handler.temp_dir)


def test_init_in_memory(memory_handler_in_memory):
    assert os.path.isdir(memory_handler_in_memory.temp_dir)
    assert memory_handler_in_memory.temp_dir.startswith(AppConfig.RAMDISK_PATH)


def test_get_temp_dir(memory_handler):
    temp_dir  = memory_handler.get_temp_dir()

    assert os.path.isdir(temp_dir)
    assert os.listdir(temp_dir) == []
    
    temp_dir_in_temp_dir = memory_handler.get_temp_dir(temp_dir)

    assert temp_dir_in_temp_dir.startswith(temp_dir)
    assert os.path.isdir(temp_dir_in_temp_dir)
    assert os.listdir(temp_dir_in_temp_dir) == []

    opened_dirs = [memory_handler.temp_dir, temp_dir, temp_dir_in_temp_dir]
    assert [type(dir) for dir in memory_handler.opened_dirs] == [TemporaryDirectory] * 3
    assert [dir.name for dir in memory_handler.opened_dirs] == opened_dirs


def test_get_temp_dir_in_memory(memory_handler_in_memory):
    temp_dir  = memory_handler_in_memory.get_temp_dir()

    assert os.path.isdir(temp_dir)
    assert os.listdir(temp_dir) == []
    assert temp_dir.startswith(AppConfig.RAMDISK_PATH)


def test_get_path(memory_handler):
    target_path = "path/to/file"

    path = memory_handler.get_path(target_path)

    assert path == target_path


def test_get_path_in_memory(memory_handler_in_memory, tmp_path):
    file_path = "path/to/file"
    source_path = os.path.join(tmp_path, file_path)
    
    os.makedirs(os.path.dirname(source_path))
    open(source_path, 'w').write("file contents")
    
    in_memory_path = memory_handler_in_memory.get_path(source_path)

    assert in_memory_path.startswith(memory_handler_in_memory.temp_dir)
    assert in_memory_path.endswith(os.path.basename(file_path))

    assert open(in_memory_path).read() == "file contents"

    open(source_path, "w").write("new file contents")

    new_in_memory_path = memory_handler_in_memory.get_path(source_path)

    assert open(new_in_memory_path).read() == "file contents"


def test_get_output_path(memory_handler):
    target_path = "path/to/file"

    path = memory_handler.get_output_path(target_path)

    assert path == target_path


def test_get_output_path_in_memory(memory_handler_in_memory, tmp_path):
    source_path_1 = os.path.join(tmp_path, "path/to/output_file_1")
    source_path_2 = os.path.join(tmp_path, "path/to/output_file_2")
    source_path_3 = os.path.join(tmp_path, "path/to/output_file_3")

    in_memory_path_1 = memory_handler_in_memory.get_output_path(source_path_1)
    in_memory_path_2 = memory_handler_in_memory.get_output_path(source_path_2, save_on_close=False)
    in_memory_path_3 = memory_handler_in_memory.get_output_path(source_path_3, save_on_close=True)

    for in_memory_path_n, source_path_n in ((in_memory_path_1, source_path_1),
                                            (in_memory_path_2, source_path_2),
                                            (in_memory_path_3, source_path_3)):
        assert os.path.isfile(in_memory_path_n)
        assert in_memory_path_n.startswith(memory_handler_in_memory.temp_dir)
        assert in_memory_path_n.endswith(os.path.basename(source_path_n))

    assert memory_handler_in_memory.saving_paths == {
        in_memory_path_1: source_path_1,
        in_memory_path_3: source_path_3
    }


def test_get_bam_path(memory_handler, tmp_path):
    bam_path = os.path.join(tmp_path, "alignment.bam")
    bai_path = os.path.join(tmp_path, "alignment.bam.bai")

    open(bam_path, "w").write("alignment contents")
    open(bai_path, "w").write("alignment index contents")

    assert memory_handler.get_bam_path(bam_path) == bam_path


def test_get_bam_path_in_memory(memory_handler_in_memory, tmp_path):
    bam_file = "alignment.bam"
    bai_file = "alignment.bam.bai"

    bam_path = os.path.join(tmp_path, bam_file)
    bai_path = os.path.join(tmp_path, bai_file)

    open(bam_path, "w").write("alignment contents")
    open(bai_path, "w").write("alignment index contents")

    in_memory_bam_path = memory_handler_in_memory.get_bam_path(bam_path)
    in_memory_bai_path = os.path.join(memory_handler_in_memory.temp_dir, bai_file)

    assert in_memory_bam_path == os.path.join(memory_handler_in_memory.temp_dir, bam_file)

    assert os.path.isfile(in_memory_bam_path)
    assert os.path.isfile(in_memory_bai_path)

    assert open(in_memory_bam_path).read() == "alignment contents"
    assert open(in_memory_bai_path).read() == "alignment index contents"


def test_save_tempfiles_to_disk(memory_handler, tmp_path):
    opened_dirs_backup = memory_handler.opened_dirs

    source_path = os.path.join(tmp_path, "output_file")
    memory_handler.get_output_path(source_path)

    memory_handler.save_tempfiles_to_disk()

    assert memory_handler.saving_paths == {}

    assert memory_handler.opened_dirs == []

    for opened_dir in opened_dirs_backup:
        assert not os.path.exists(opened_dir.name)


def test_save_tempfiles_to_disk_in_memory(memory_handler_in_memory, tmp_path):
    opened_dirs_backup = memory_handler_in_memory.opened_dirs

    source_path_1 = os.path.join(tmp_path, "output_file_1")
    source_path_2 = os.path.join(tmp_path, "output_file_2")
    source_path_3 = os.path.join(tmp_path, "output_file_3")

    in_memory_path_1 = memory_handler_in_memory.get_output_path(source_path_1)
    in_memory_path_2 = memory_handler_in_memory.get_output_path(source_path_2, save_on_close=False)
    in_memory_path_3 = memory_handler_in_memory.get_output_path(source_path_3, save_on_close=True)

    open(in_memory_path_1, "w").write("file contents 1")
    open(in_memory_path_2, "w").write("file contents 2")
    open(in_memory_path_3, "w").write("file contents 3")

    memory_handler_in_memory.save_tempfiles_to_disk()

    assert memory_handler_in_memory.saving_paths == {}

    assert open(source_path_1).read() == "file contents 1"
    assert open(source_path_3).read() == "file contents 3"

    assert not os.path.exists(source_path_2)

    assert memory_handler_in_memory.opened_dirs == []

    for opened_dir in opened_dirs_backup:
        assert not os.path.exists(opened_dir.name)


def test_close_without_saving(memory_handler, tmp_path):
    opened_dirs_backup = memory_handler.opened_dirs

    source_path = os.path.join(tmp_path, "output_file")
    memory_handler.get_output_path(source_path)

    memory_handler.close_without_saving()

    assert memory_handler.saving_paths == {}

    assert memory_handler.opened_dirs == []

    for opened_dir in opened_dirs_backup:
        assert not os.path.exists(opened_dir.name)


def test_close_without_saving_in_memory(memory_handler_in_memory, tmp_path):
    opened_dirs_backup = memory_handler_in_memory.opened_dirs

    source_path_1 = os.path.join(tmp_path, "output_file_1")
    source_path_2 = os.path.join(tmp_path, "output_file_2")
    source_path_3 = os.path.join(tmp_path, "output_file_3")

    memory_handler_in_memory.get_output_path(source_path_1)
    memory_handler_in_memory.get_output_path(source_path_2, save_on_close=False)
    memory_handler_in_memory.get_output_path(source_path_3, save_on_close=True)

    memory_handler_in_memory.close_without_saving()

    assert memory_handler_in_memory.saving_paths == {}

    for source_path_n in (source_path_1, source_path_2, source_path_3):
        assert not os.path.exists(source_path_n)

    assert memory_handler_in_memory.opened_dirs == []

    for opened_dir in opened_dirs_backup:
        assert not os.path.exists(opened_dir.name)


def test_context_management(memory_handler, mocker):
    with memory_handler:
        assert os.path.isdir(memory_handler.temp_dir)
        assert len(memory_handler.opened_dirs) == 1
        assert memory_handler.opened_dirs[0].name == memory_handler.temp_dir

    assert memory_handler.temp_dir == None
    assert memory_handler.opened_dirs == []

    with memory_handler:
        assert os.path.isdir(memory_handler.temp_dir)
        assert len(memory_handler.opened_dirs) == 1
        assert memory_handler.opened_dirs[0].name == memory_handler.temp_dir

    assert memory_handler.temp_dir == None
    assert memory_handler.opened_dirs == []

    with pytest.raises(Exception):
        with memory_handler:
            assert os.path.isdir(memory_handler.temp_dir)
            assert len(memory_handler.opened_dirs) == 1
            assert memory_handler.opened_dirs[0].name == memory_handler.temp_dir
            raise Exception("Test exception")

    assert memory_handler.temp_dir == None
    assert memory_handler.opened_dirs == []


def test_context_management_in_memory(memory_handler_in_memory, mocker):
    mocker.patch.object(memory_handler_in_memory, "get_temp_dir")
    mocker.patch.object(memory_handler_in_memory, "save_tempfiles_to_disk")
    mocker.patch.object(memory_handler_in_memory, "close_without_saving")

    with memory_handler_in_memory:
        pass
    
    assert memory_handler_in_memory.save_tempfiles_to_disk.call_count == 1

    with memory_handler_in_memory:
        pass
    
    assert memory_handler_in_memory.save_tempfiles_to_disk.call_count == 2

    with pytest.raises(Exception):
        with memory_handler_in_memory:
            raise Exception("Test exception")
        
    assert memory_handler_in_memory.close_without_saving.call_count == 1
