
import os

from IPython import embed

import numpy

from matplotlib import pyplot

from mangadap.contrib import xjmc

def test_component_select():
    comp = numpy.array([ 0,  1,  2,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 14,
                        15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,  8])

    valid = numpy.ones(comp.size, dtype=bool)

    _comp, comp_map = xjmc._reset_components(comp, valid)
    assert numpy.array_equal(comp, _comp), 'Component list should not have changed.'

    # Remove one template that doesn't cause the removal of an entire
    # component.
    valid[3] = False
    _comp, comp_map = xjmc._reset_components(comp, valid)
    assert numpy.array_equal(comp[valid], _comp), 'Component IDs should not have changed'
    valid[3] = True

    # Remove a template that constitutes an entire component
    valid[12] = False
    _comp, comp_map = xjmc._reset_components(comp, valid)
    assert len(_comp) == len(comp)-1, 'Component was not removed'
    assert numpy.array_equal(_comp[12:-1]+1, comp[13:-1]), \
            'Components were not reordered correctly.'
    valid[12] = True


def test_constraint_select():
    comp = numpy.array([ 0, 1, 2, 2, 3, 4])
    valid = numpy.ones(comp.size, dtype=bool)

    A_ineq = numpy.array([[ 0., 0.91, 0., 0., 0., 0., 0., -1., 0., 0.],
                          [ 0., -1.1, 0., 0., 0., 0., 0.,  1., 0., 0.]])
    b_ineq = numpy.zeros(2, dtype=float)

    # Removes a template that has other templates in this component, so the
    # return arrays should be identical.
    valid[2] = False
    _A, _b = xjmc._reset_kinematic_constraints(comp, valid, A_ineq, b_ineq)
    assert numpy.array_equal(A_ineq, _A) and numpy.array_equal(b_ineq, _b), \
            'No components should have been removed'
    valid[2] = True

    # Removes the only template in a component, but the component is not
    # constrained by A,b. So the shape of the return A matrix should have the
    # relevant component columns removed.
    valid[1] = False
    _A, _b = xjmc._reset_kinematic_constraints(comp, valid, A_ineq, b_ineq)
    assert numpy.array_equal(A_ineq[:,numpy.repeat([True, False, True, True, True],2)], _A) \
                and numpy.array_equal(b_ineq, _b), 'No components should have been removed'
    valid[1] = True

    # Removes a template constituting an entire component constrained by A,b,
    # such that no constraints remain
    valid[0] = False
    _A, _b = xjmc._reset_kinematic_constraints(comp, valid, A_ineq, b_ineq)
    assert _A is None and _b is None, 'No constraints should remain'
    valid[0] = True

    # Try with more constraints    
    A_ineq = numpy.array([[ 0., 0.91, 0., 0., 0.,   0., 0., -1., 0., 0.],
                          [ 0., -1.1, 0., 0., 0.,   0., 0.,  1., 0., 0.],
                          [ 0.,   0., 0., 0., 0., 0.91, 0., -1., 0., 0.],
                          [ 0.,   0., 0., 0., 0., -1.1, 0.,  1., 0., 0.]])
    b_ineq = numpy.zeros(4, dtype=float)

    # Removes the only template in a component, but the component is not
    # constrained by A,b. So the shape of the return A matrix should have the
    # relevant component columns removed.
    valid[1] = False
    _A, _b = xjmc._reset_kinematic_constraints(comp, valid, A_ineq, b_ineq)
    assert numpy.array_equal(A_ineq[:,numpy.repeat([True, False, True, True, True],2)], _A) \
                and numpy.array_equal(b_ineq, _b), 'No components should have been removed'
    valid[1] = True

    # Removes a template constituting an entire component constrained by 2 of
    # the 4 constraints defined by A,b, such only 2 remain in the output.
    valid[0] = False
    _A, _b = xjmc._reset_kinematic_constraints(comp, valid, A_ineq, b_ineq)
    assert numpy.array_equal(A_ineq[2:,numpy.repeat([False, True, True, True, True],2)], _A) \
                and numpy.array_equal(b_ineq[2:], _b), 'No components should have been removed'
    valid[0] = True

    # Removes a template constituting an entire component that included in all
    # the constraints defined by A,b, such that no constraints remain
    valid[4] = False
    _A, _b = xjmc._reset_kinematic_constraints(comp, valid, A_ineq, b_ineq)
    assert _A is None and _b is None, 'No constraints should remain'
    valid[4] = True


