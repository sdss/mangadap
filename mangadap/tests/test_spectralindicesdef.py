
import pytest

from mangadap.proc.spectralindices import available_spectral_index_databases
from mangadap.par.artifactdb import ArtifactDB
from mangadap.par.absorptionindexdb import AbsorptionIndexDB
from mangadap.par.bandheadindexdb import BandheadIndexDB

def test_avail():
    database_list = available_spectral_index_databases()
    assert len(database_list) > 0, 'No spectral index databases available'

    # For the available databases, make sure that the ancillary
    # databases can be loaded.
    for database in database_list:
        if database['artifacts'] is not None:
            artdb = ArtifactDB.from_key(database['artifacts'])
        if database['absindex'] is not None:
            absdb = AbsorptionIndexDB.from_key(database['absindex'])
        if database['bandhead'] is not None:
            bhddb = BandheadIndexDB.from_key(database['bandhead'])

