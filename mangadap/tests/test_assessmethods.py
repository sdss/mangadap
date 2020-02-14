
import pytest

from mangadap.proc.reductionassessments import available_reduction_assessments

def test_avail():
    method_list = available_reduction_assessments()
    assert len(method_list) > 0, 'No reduction assessment methods available'

