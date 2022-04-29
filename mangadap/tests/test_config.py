from pathlib import Path

from IPython import embed

from mangadap.config.manga import MaNGAConfig
from mangadap.tests.util import data_test_file, remote_data_file, requires_remote


def test_init():
    cfg = MaNGAConfig(7815, 3702)
    assert isinstance(cfg.directory_path, Path), 'Directory should be a Path instance'
    assert cfg.plate == 7815, 'Plate changed'
    assert cfg.log, 'Log binning should be true'
    assert cfg.mode == 'CUBE', 'Should be in CUBE mode'
    assert cfg.file_name == 'manga-7815-3702-LOGCUBE.fits.gz'


def test_copy():
    cfg = MaNGAConfig(7815, 3702)
    _cfg = cfg.copy()
    assert cfg.file_name == _cfg.file_name, 'Bad name copy'
    _cfg.mode = 'RSS'
    assert cfg.file_name != _cfg.file_name, 'Bad name change'


def test_from_config():
    cfg = MaNGAConfig.from_config(data_test_file("datacube.ini"))
    assert isinstance(cfg.directory_path, Path), 'Directory should be a Path instance'
    # TODO: This passes locally, but will fail the tox CI tests.  Determine how
    # to perform this in tox?
    #assert cfg.directory_path == Path(remote_data_file()).resolve(), 'Directory path changed'
    assert cfg.plate == 7815, 'Plate changed'
    assert cfg.log, 'Log binning should be true'
    assert cfg.mode == 'CUBE', 'Should be in CUBE mode'
    assert cfg.file_name == 'manga-7815-3702-LOGCUBE.fits.gz'


@requires_remote
def test_from_file():
    file = remote_data_file('manga-7815-3702-LOGCUBE.fits.gz')
    cfg = MaNGAConfig.from_file(file)

    assert isinstance(cfg.directory_path, Path), 'Directory should be a Path instance'
    assert cfg.directory_path == Path(remote_data_file()).resolve(), 'Directory path changed'
    assert cfg.plate == 7815, 'Plate changed'
    assert cfg.log, 'Log binning should be true'
    assert cfg.mode == 'CUBE', 'Should be in CUBE mode'
    assert cfg.file_name == 'manga-7815-3702-LOGCUBE.fits.gz'


