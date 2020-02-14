
import pytest
import os
from IPython import embed

import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

from mangadap.proc.emissionlinemoments import available_emission_line_moment_databases

from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.emissionmomentsdb import EmissionMomentsDB

def test_avail():
    # Get the available methods
    db_list = available_emission_line_moment_databases()
    assert len(db_list) > 0, 'No emission-line-moment databases available'

    # For the available methods, make sure that the ancillary databases
    # and templates can be loaded.
    for db in db_list:
        print(db['key'])
        if db['artifacts'] is not None:
            artdb = ArtifactDB.from_key(db['artifacts'])
        if db['passbands'] is not None:
            emomdb = EmissionMomentsDB.from_key(db['passbands'])

