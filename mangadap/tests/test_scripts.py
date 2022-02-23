
import pytest
import os
import subprocess
import shutil

from IPython import embed

from mangadap.tests.util import remote_data_file, data_test_file, drp_test_version
from mangadap.tests.util import requires_remote, requires_drpcomplete, requires_drpall

from mangadap.scripts.calculate_covariance import CalculateCovariance
from mangadap.scripts.manga_dap import MangaDap
from mangadap.scripts.write_dap_config import WriteDapConfig
from mangadap.util.parser import DefaultConfig


def script_help_okay(executable):
    # Basic help call
    return subprocess.call([executable, '-h'], stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL) == 0

@requires_remote
def test_calculate_covariance():

    assert script_help_okay('dap_calculate_covariance'), 'Basic help call failed'

    ofile = 'test.fits'
    if os.path.isfile(ofile):
        # Clean-up a previous failure
        os.remove(ofile)

    # Defaults to central channel
    # TODO: Test contents?
    CalculateCovariance.main(CalculateCovariance.parse_args(['7815', '3702', ofile, '-d',
                                                             remote_data_file()]))
    assert os.path.isfile(ofile), 'Output file not written.'

    # Do a number of channels
    os.remove(ofile)
    CalculateCovariance.main(CalculateCovariance.parse_args(['7815', '3702', ofile, '-n', '11',
                                                             '-d', remote_data_file()]))
    assert os.path.isfile(ofile), 'Output file not written.'

    # Run a specific wavelength
    os.remove(ofile)
    CalculateCovariance.main(CalculateCovariance.parse_args(['7815', '3702', ofile, '-w',
                                                             '4500.', '-d', remote_data_file()]))
    assert os.path.isfile(ofile), 'Output file not written.'

    # Clean-up
    os.remove(ofile)


# TODO: Add some remote files to test this?
def test_construct_dapall():
    assert script_help_okay('construct_dapall'), 'Basic help call failed'


# TODO: Include this in the construct_dapall test and test the
# constructed file?
def test_dapall_qa():
    assert script_help_okay('dapall_qa'), 'Basic help call failed'


# TODO: Include a DRPComplete file in remote?
def test_find_repeat_observations():
    assert script_help_okay('dap_find_repeat_observations'), 'Basic help call failed'


# TODO: Add some remote files?
def test_fit_residuals():
    assert script_help_okay('dap_fit_residuals'), 'Basic help call failed'


# TODO: Add some remote files?
def test_inspector():
    assert script_help_okay('manga_dap_inspector'), 'Basic help call failed'


def test_manga_dap_import():
    with pytest.raises(ImportError):
        MangaDap.main(MangaDap.parse_args(['-c', data_test_file('datacube.ini'), '-m', 'junk']))
    
    with pytest.raises(AttributeError):
        MangaDap.main(MangaDap.parse_args(['-c', data_test_file('datacube.ini'), '-o', 'junk']))
    

@requires_remote
def test_manga_dap():
    odir = data_test_file('manga_dap_output')

    # Clean up previous failure
    if os.path.isdir(odir):
        shutil.rmtree(odir)

    # Run the DAP. The binning in plan.par is set to ALL binning, so
    # this run of the DAP just analyzes one spectrum and takes about a
    # minute.
    MangaDap.main(MangaDap.parse_args(['-c', data_test_file('datacube.ini'),
                                       '-p', data_test_file('plan.par'), '-a', odir]))

    # Re-run to use existing files.  Takes about 40s.
    MangaDap.main(MangaDap.parse_args(['-c', data_test_file('datacube.ini'),
                                       '-p', data_test_file('plan.par'), '-a', odir]))

    # Clean up
    shutil.rmtree(odir)


# TODO: Add some remote data?
def test_plate_fit_qa():
    assert script_help_okay('dap_plate_fit_qa'), 'Basic help call failed'


# TODO: Add some remote data?
def test_ppxffit_qa():
    assert script_help_okay('dap_ppxffit_qa'), 'Basic help call failed'


# TODO: Add some remote data?
def test_spotcheck():
    assert script_help_okay('spotcheck_dap_maps'), 'Basic help call failed'


# TODO: Add some remote data?
def test_template_flux_norm():
    assert script_help_okay('dap_template_flux_norm'), 'Basic help call failed'


# TODO: Add a remote DRPComplete file!
def test_write_dap_config_basic():
    assert script_help_okay('write_dap_config'), 'Basic help call failed'


@requires_drpcomplete
def test_write_dap_config_drpcomplete():
    ofile = 'test.ini'
    if os.path.isfile(ofile):
        # Clean-up a previous failure
        os.remove(ofile)

    drpc = remote_data_file('drpcomplete_{0}.fits'.format(drp_test_version))
    WriteDapConfig.main(WriteDapConfig.parse_args(['7815', '3702', ofile, '-c', drpc]))

    assert os.path.isfile(ofile), 'Output file not written.'

    cfg = DefaultConfig(ofile)
    assert cfg.getint('plate') == 7815, 'Wrong plate number'
    assert cfg.get('directory_path') is None, 'No directory path should be defined.'
    assert cfg.getfloat('z') == 0.0293823, 'Bad redshift'

    # Clean-up a previous failure
    os.remove(ofile)


@requires_drpall
def test_write_dap_config_drpall():
    ofile = 'test.ini'
    if os.path.isfile(ofile):
        # Clean-up a previous failure
        os.remove(ofile)

    drpall = remote_data_file('drpall-{0}.fits'.format(drp_test_version))
    WriteDapConfig.main(WriteDapConfig.parse_args(['7815', '3702', ofile, '-a', drpall]))

    assert os.path.isfile(ofile), 'Output file not written.'

    cfg = DefaultConfig(ofile)
    assert cfg.getint('plate') == 7815, 'Wrong plate number'
    assert cfg.get('directory_path') is None, 'No directory path should be defined.'
    assert cfg.getfloat('z') == 0.0293823, 'Bad redshift'

    # Clean-up a previous failure
    os.remove(ofile)


# TODO: Add some remote files?
def test_status():
    assert script_help_okay('dap_status'), 'Basic help call failed'


# TODO: Improve this!
def test_run():
    assert script_help_okay('rundap'), 'Basic help call failed'

