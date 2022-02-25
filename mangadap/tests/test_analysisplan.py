from IPython import embed

from mangadap.tests.util import data_test_file
from mangadap.par.analysisplan import AnalysisPlanSet


def test_default():
    plan = AnalysisPlanSet.default()
    assert len(plan) == 1, 'Default is to perform one plan.'
    assert plan[0]['drpqa_key'] == 'SNRG', 'Default DRP reduction QA key changed.'
    assert plan[0]['elfit_key'] == 'EFITMPL11HCDB', 'Default emission-line fit key changed.'


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

