from pathlib import Path
from IPython import embed

from mangadap.tests.util import data_test_file
from mangadap.config.analysisplan import AnalysisPlan
from mangadap.config.manga import MaNGAAnalysisPlan, MaNGAConfig


def test_default():
    plan = AnalysisPlan.default()
    assert len(plan.keys()) == 1, 'Default is to perform one plan.'
    assert list(plan.keys())[0] == 'default', 'Name changed'

    assert plan.rdxqa['default']['key'] == 'SNRG', 'Default DRP reduction QA key changed.'
    assert plan.elfit['default']['key'] == 'EFITHC2', 'Default emission-line fit key changed.'
    assert plan.common_path() == Path('.').resolve() / 'common', 'Common path changed'


def test_read():
    plan = AnalysisPlan.from_toml(data_test_file('global_bin.toml'))

    assert len(plan.keys()) == 1, 'Number of example plans changed'
    assert list(plan.keys())[0] == 'default', 'Name changed'

    assert plan.binning['default']['key'] == 'ALL', 'Default DRP reduction QA key changed.'
    assert plan.elfit['default']['key'] == 'EFITSSP', 'Default emission-line fit key changed.'


def test_multi_read():
    plan = AnalysisPlan.from_toml(data_test_file('dr17.toml'))

    assert len(plan.keys()) == 4, 'Number of example plans changed'

    assert plan['plan3']['key'] == 'HYB10-MILESHC-MASTARSSP', 'Plan key changed'
    assert plan.binning['plan1']['key'] == 'SPX', 'Binning changed'
    assert plan.elfit['plan4']['key'] == 'EFITHC2DB', 'Emission-line fitting key changed'


def test_mangaplan():
    cfg = MaNGAConfig(7815, 3702)
    plan = MaNGAAnalysisPlan.from_toml(data_test_file('dr17.toml'), cube=cfg)

    assert plan.dap_file_root(cfg) == 'manga-7815-3702', 'Bad root'
    assert plan.dap_file_root(cfg, mode='MAPS', plan_index=0) \
                == 'manga-7815-3702-MAPS-SPX-MILESHC-MASTARSSP', 'Bad full root'
    assert plan.common_path().parts[-2:] == (str(cfg.plate), str(cfg.ifudesign)), \
                'Bad common subdirectories'
    assert plan.method_path().parts[-2:] == (str(cfg.plate), str(cfg.ifudesign)), \
                'Bad method subdirectories'
    assert plan.method_path(qa=True).parts[-2:] == (str(cfg.ifudesign), 'qa'), \
                'Bad qa subdirectories'


