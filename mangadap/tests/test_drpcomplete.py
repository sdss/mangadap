import pytest

import os

from IPython import embed

from mangadap.survey.drpcomplete import DRPComplete
from mangadap.util.parser import DefaultConfig
from mangadap.tests.util import remote_data_file, drp_test_version, requires_drpcomplete

@requires_drpcomplete
def test_write_cfg():
    ofile = 'test.ini'
    if os.path.isfile(ofile):
        # Remove existing files from failed tests
        os.remove(ofile)

    # Read the test DRPComplete file.
    # TODO: Temporarily override
    drp_test_version = 'v3_1_1'
    drpc = DRPComplete(drpver=drp_test_version, directory_path=remote_data_file(), readonly=True)

    # Write the base-level configuration file
    drpc.write_config(ofile, plate=7815, ifudesign=3702)

    # Read it and check the output
    cfg = DefaultConfig(ofile)
    assert cfg.getint('plate') == 7815, 'Plate number is wrong'
    assert cfg.getbool('log'), 'Should be selecting the logarithmically binned data'
    assert cfg.get('sres_ext') is None, 'Spectral resolution extension should be undefined'
    assert cfg.getfloat('z') == 2.9382300e-02, 'Bad redshift'

    # Try to write it again with overwrite set to False
    with pytest.raises(FileExistsError):
        drpc.write_config(ofile, plate=7815, ifudesign=3702, overwrite=False)

    # Set the spectral resolution flags
    drpc.write_config(ofile, plate=7815, ifudesign=3702, overwrite=True, sres_ext='SPECRES',
                      sres_fill=False)

    # Read it and check the output
    cfg = DefaultConfig(ofile)
    assert cfg.getint('plate') == 7815, 'Plate number is wrong'
    assert cfg.getbool('log'), 'Should be selecting the logarithmically binned data'
    assert cfg.get('sres_ext') == 'SPECRES', 'Spectral resolution extension incorrect'
    assert cfg.getfloat('z') == 2.9382300e-02, 'Bad redshift'

    # Clean-up
    os.remove(ofile)


