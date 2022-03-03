from pathlib import Path
from IPython import embed

from mangadap.tests.util import data_test_file
from mangadap.config.analysisplan import AnalysisPlanSet
from mangadap.config.manga import MaNGAAnalysisPlan, MaNGAConfig


def test_default():
    plan = AnalysisPlanSet.default()
    assert len(plan) == 1, 'Default is to perform one plan.'
    assert plan[0]['drpqa_key'] == 'SNRG', 'Default DRP reduction QA key changed.'
    assert plan[0]['elfit_key'] == 'EFITMPL11HCDB', 'Default emission-line fit key changed.'
    assert plan.common_path() == Path('.').resolve() / 'common', 'Common path changed'


def test_read():
    plan = AnalysisPlanSet.from_par_file(data_test_file('plan.par'))

    assert len(plan) == 1, 'Number of example plans changed'
    assert plan[0]['bin_key'] == 'ALL', 'Binning changed'
    assert plan[0]['elfit_key'] == 'EFITMPL11HC', 'Emission-line fitting key changed'


def test_multi_read():
    plan = AnalysisPlanSet.from_par_file(data_test_file('dr17_plan.par'))

    assert len(plan) == 4, 'Number of example plans changed'
    assert plan[0]['bin_key'] == 'SPX', 'Binning changed'
    assert plan[-1]['elfit_key'] == 'EFITMPL11HCDB', 'Emission-line fitting key changed'
    assert plan[2]['key'] == 'HYB10-MILESHC-MASTARSSP', 'Plan key changed'


def test_mangaplan():
    cfg = MaNGAConfig(7815, 3702)
    plan = MaNGAAnalysisPlan.from_par_file(data_test_file('dr17_plan.par'), cube=cfg)

    assert plan.dap_file_root(cfg) == 'manga-7815-3702', 'Bad root'
    assert plan.dap_file_root(cfg, mode='MAPS', plan_index=0) \
                == 'manga-7815-3702-MAPS-SPX-MILESHC-MASTARSSP', 'Bad full root'
    assert plan.common_path().parts[-2:] == (str(cfg.plate), str(cfg.ifudesign)), \
                'Bad common subdirectories'
    assert plan.method_path().parts[-2:] == (str(cfg.plate), str(cfg.ifudesign)), \
                'Bad method subdirectories'
    assert plan.method_path(qa=True).parts[-2:] == (str(cfg.ifudesign), 'qa'), \
                'Bad qa subdirectories'



