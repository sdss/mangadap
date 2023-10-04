from IPython import embed

import numpy

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



def test_component_setup_no_stellar():

    component = numpy.array([0, 0, 1, 2, 2])
    moments = numpy.array([2, 2, 2])
    npar = numpy.sum(moments)
    gas_template = numpy.array([True, True, True, True, True])
    start = numpy.zeros((2,npar), dtype=float)
    start[0] = numpy.concatenate([[ 0., 100.],
                                  [10., 200.],
                                  [20., 300.]])
    start[1] = numpy.concatenate([[50., 150.],
                                  [60., 250.],
                                  [70., 350.]])
    fixed = numpy.concatenate([[False, False],
                               [False, False],
                               [False, False]])
    gas_start = numpy.array([[30., 100.],
                             [80., 150.]])

    # Basic return
    _component, _moments, _start, _fixed \
        = xjmc._ppxf_component_setup(component, moments, gas_template, start, fixed)
    assert _component is not component, 'Should be a copy'
    assert numpy.all(_component == component), 'Should be identical'
    
    # Set single gas component
    _component, _moments, _start, _fixed \
        = xjmc._ppxf_component_setup(component, moments, gas_template, start, fixed,
                                     single_gas_component=True)
    assert numpy.all(_component == 0), 'Should be one component'
    assert _moments.size == 1, 'Should only be one moment'
    assert _moments[0] == 2, 'Moments should be Gaussian'
    assert _start[0].tolist() == [10., 200.], 'Bad mean start'
    assert not numpy.any(_fixed), 'Nothing should be fixed'

    # Set single gas component with provided start
    _component, _moments, _start, _fixed \
        = xjmc._ppxf_component_setup(component, moments, gas_template, start, fixed,
                                     single_gas_component=True, gas_start=gas_start)
    assert numpy.all(_component == 0), 'Should be one component'
    assert _moments.size == 1, 'Should only be one moment'
    assert _moments[0] == 2, 'Moments should be Gaussian'
    assert numpy.all(_start == gas_start), 'Bad start'

    # Provide start for all gas components
    _component, _moments, _start, _fixed \
        = xjmc._ppxf_component_setup(component, moments, gas_template, start, fixed,
                                     gas_start=gas_start)
    assert _component is not component, 'Should be a copy'
    assert numpy.all(_component == component), 'Should be identical'
    assert numpy.all(_start == numpy.tile(gas_start, (3,))), 'All gas starts should be the same'


def test_component_setup_one_stellar():

    component = numpy.array([0, 0, 1, 2, 2])
    # Gaussian stellar moments
    moments = numpy.array([2, 2, 2])
    npar = numpy.sum(moments)
    gas_template = numpy.array([False, False, True, True, True])
    start = numpy.zeros((2,npar), dtype=float)
    start[0] = numpy.concatenate([[ 0., 100.],
                                  [10., 200.],
                                  [20., 300.]])
    start[1] = numpy.concatenate([[50., 150.],
                                  [60., 250.],
                                  [70., 350.]])
    fixed = numpy.concatenate([[True, True],
                               [False, False],
                               [False, False]])
    gas_start = numpy.array([[30., 100.],
                             [80., 150.]])

    # Basic return
    _component, _moments, _start, _fixed \
        = xjmc._ppxf_component_setup(component, moments, gas_template, start, fixed)
    assert _component is not component, 'Should be a copy'
    assert numpy.all(_component == component), 'Should be identical'
    
    # Set single gas component
    _component, _moments, _start, _fixed \
        = xjmc._ppxf_component_setup(component, moments, gas_template, start, fixed,
                                     single_gas_component=True)
    assert numpy.amax(_component) == 1, 'Should only be two components 0,1'
    assert len(_moments) == 2, 'Should only be two components'
    assert numpy.all(start[:,:2] == _start[:,:2]), 'Should copy stellar start'
    assert numpy.all(_start[:,2:] == numpy.array([[15., 250.],
                                                  [65., 300.]])), 'Bad mean start'

    # Set single gas component with provided start
    _component, _moments, _start, _fixed \
        = xjmc._ppxf_component_setup(component, moments, gas_template, start, fixed,
                                     single_gas_component=True, gas_start=gas_start)
    assert numpy.amax(_component) == 1, 'Should only be two components 0,1'
    assert len(_moments) == 2, 'Should only be two components'
    assert numpy.all(start[:,:2] == _start[:,:2]), 'Should copy stellar start'
    assert numpy.all(_start[:,2:] == gas_start), 'Bad start copy'

    # Provide start for all gas components
    _component, _moments, _start, _fixed \
        = xjmc._ppxf_component_setup(component, moments, gas_template, start, fixed,
                                     gas_start=gas_start)
    assert _component is not component, 'Should be a copy'
    assert numpy.all(_component == component), 'Should be identical'
    assert numpy.all(_start[:,2:] == numpy.tile(gas_start, (2,))), 'All gas starts should be the same'

    # 4 stellar moments
    moments = numpy.array([4, 2, 2])
    npar = numpy.sum(moments)
    start = numpy.zeros((2,npar), dtype=float)
    start[0] = numpy.concatenate([[ 0., 100., 0.1, 0.1],
                                  [10., 200.],
                                  [20., 300.]])
    start[1] = numpy.concatenate([[50., 150., 0.1, 0.1],
                                  [60., 250.],
                                  [70., 350.]])
    fixed = numpy.concatenate([[True, True, True, True],
                               [False, False],
                               [False, False]])

    # Basic return
    _component, _moments, _start, _fixed \
        = xjmc._ppxf_component_setup(component, moments, gas_template, start, fixed)
    assert _component is not component, 'Should be a copy'
    assert numpy.all(_component == component), 'Should be identical'
    
    # Set single gas component
    _component, _moments, _start, _fixed \
        = xjmc._ppxf_component_setup(component, moments, gas_template, start, fixed,
                                     single_gas_component=True)
    assert numpy.amax(_component) == 1, 'Should only be two components 0,1'
    assert len(_moments) == 2, 'Should only be two components'
    assert numpy.all(start[:,:4] == _start[:,:4]), 'Should copy stellar start'
    assert numpy.all(_start[:,4:] == numpy.array([[15., 250.],
                                                  [65., 300.]])), 'Bad mean start'


def test_component_setup_two_stellar():

    component = numpy.array([0, 0, 1, 2, 2, 3])
    # Gaussian stellar moments
    moments = numpy.array([2, 2, 2, 2])
    npar = numpy.sum(moments)
    gas_template = numpy.array([False, False, False, True, True, True])
    start = numpy.zeros((2,npar), dtype=float)
    start[0] = numpy.concatenate([[ 0., 100.],
                                  [10., 200.],
                                  [20., 300.],
                                  [30., 400.]])
    start[1] = numpy.concatenate([[50., 150.],
                                  [60., 250.],
                                  [70., 350.],
                                  [80., 450.]])
    fixed = numpy.concatenate([[True, True],
                               [True, True],
                               [False, False],
                               [False, False]])
    gas_start = numpy.array([[30., 100.],
                             [80., 150.]])

    # Basic return
    _component, _moments, _start, _fixed \
        = xjmc._ppxf_component_setup(component, moments, gas_template, start, fixed)
    assert _component is not component, 'Should be a copy'
    assert numpy.all(_component == component), 'Should be identical'

    # Set single gas component
    _component, _moments, _start, _fixed \
        = xjmc._ppxf_component_setup(component, moments, gas_template, start, fixed,
                                     single_gas_component=True)
    assert numpy.amax(_component) == 2, 'Should be three components 0,1,2'
    assert len(_moments) == 3, 'Should be three components'
    assert numpy.all(start[:,:4] == _start[:,:4]), 'Should copy stellar start'
    assert numpy.all(_start[:,4:] == numpy.array([[25., 350.],
                                                  [75., 400.]])), 'Bad mean start'


def test_reset_components():

    c = numpy.arange(5)
    valid = numpy.array([0, 1, 1, 0, 1], dtype=bool)
    _c, c_map = xjmc._reset_components(c, valid)

    assert numpy.array_equal(_c, numpy.arange(numpy.sum(valid))), 'Bad renumbering'
    assert numpy.array_equal(c[valid], c_map), 'Bad mapping'







