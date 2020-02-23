
import pytest
import os
import subprocess


from IPython import embed

from mangadap.tests.util import remote_data_file, data_test_file, drp_test_version
from mangadap.tests.util import requires_remote, requires_drpcomplete, requires_drpall

from mangadap.scripts import calculate_covariance
from mangadap.scripts import manga_dap
from mangadap.scripts import write_dap_config
from mangadap.util.parser import DefaultConfig

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

def script_help_okay(executable):
    # Basic help call
    logfile = 'test.log'
    if os.path.isfile(logfile):
        os.remove(logfile)
    with open(logfile, 'w') as f:
        retval = subprocess.call([executable, '-h'], stdout=f, stderr=f)
    os.remove(logfile)
    return retval == 0


@requires_remote
def test_calculate_covariance():

    assert script_help_okay('dap_calculate_covariance'), 'Basic help call failed'

    ofile = 'test.fits'
    if os.path.isfile(ofile):
        # Clean-up a previous failure
        os.remove(ofile)
    # Defaults to central channel
    # TODO: Test contents?
    calculate_covariance.main(calculate_covariance.parse_args(['7815', '3702', ofile, '-d',
                                                               remote_data_file()]))
    assert os.path.isfile(ofile), 'Output file not written.'

    # Do a number of channels
    os.remove(ofile)
    calculate_covariance.main(calculate_covariance.parse_args(['7815', '3702', ofile, '-n', '11',
                                                               '-d', remote_data_file()]))
    assert os.path.isfile(ofile), 'Output file not written.'

    # Run a specific wavelength
    os.remove(ofile)
    calculate_covariance.main(calculate_covariance.parse_args(['7815', '3702', ofile, '-w',
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


# TODO: Sort out an execution mode that can run in a few minutes.
def test_manga_dap():
    assert script_help_okay('manga_dap'), 'Basic help call failed'


def test_manga_dap():
    assert script_help_okay('manga_dap'), 'Basic help call failed'


def test_manga_dap_import():
    with pytest.raises(ImportError):
        manga_dap.main(manga_dap.parse_args(['-c', data_test_file('datacube.ini'), '-m', 'junk']))
    
    with pytest.raises(AttributeError):
        manga_dap.main(manga_dap.parse_args(['-c', data_test_file('datacube.ini'), '-o', 'junk']))
    

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
def test_write_dap_par():
    assert script_help_okay('write_dap_par'), 'Basic help call failed'


def test_write_dap_config_basic():
    assert script_help_okay('write_dap_config'), 'Basic help call failed'


@requires_drpcomplete
def test_write_dap_config_drpcomplete():
    ofile = 'test.ini'
    if os.path.isfile(ofile):
        # Clean-up a previous failure
        os.remove(ofile)

    drpc = remote_data_file('drpcomplete_{0}.fits'.format(drp_test_version))
    write_dap_config.main(write_dap_config.parse_args(['7815', '3702', ofile, '-c', drpc]))

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
    write_dap_config.main(write_dap_config.parse_args(['7815', '3702', ofile, '-a', drpall]))

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


if __name__ == '__main__':
    test_manga_dap_import()



