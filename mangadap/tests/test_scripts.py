
import pytest
import os
import subprocess

from IPython import embed

from mangadap.tests.util import remote_data_file, requires_remote

from mangadap.scripts import calculate_covariance

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

@requires_remote
def test_calculate_covariance():

    # Basic help call
    logfile = 'test.log'
    if os.path.isfile(logfile):
        os.remove(logfile)
    with open(logfile, 'w') as f:
        retval = subprocess.call(['dap_calculate_covariance', '-h'], stdout=f, stderr=f)
    os.remove(logfile)
    assert retval == 0, 'Basic help call failed'

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
    # Basic help call
    logfile = 'test.log'
    if os.path.isfile(logfile):
        os.remove(logfile)
    with open(logfile, 'w') as f:
        retval = subprocess.call(['construct_dapall', '-h'], stdout=f, stderr=f)
    os.remove(logfile)
    assert retval == 0, 'Basic help call failed'

# TODO: Include this in the construct_dapall test and test the
# constructed file?
def test_dapall_qa():
    # Basic help call
    logfile = 'test.log'
    if os.path.isfile(logfile):
        os.remove(logfile)
    with open(logfile, 'w') as f:
        retval = subprocess.call(['dapall_qa', '-h'], stdout=f, stderr=f)
    os.remove(logfile)
    assert retval == 0, 'Basic help call failed'

if __name__ == '__main__':
    test_dapall_qa()

