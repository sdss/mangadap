
Construction of the ``DAPTYPE`` is based on the keyword strings for
the ``binning_method``, the ``stellar_continuum_templates``, and the
``emission_line_model_templates``, with the strings separated by
dashes. For example, the ``DAPTYPE`` would be
``'HYB10-MILESHC-MILESHC'`` when the binning method is ``HYB10`` and
the ``MILESHC`` library is used for the continuum templates in both
the stellar and emission-line fitting modules.

You can construct a list of all the ``DAPTYPE`` strings for the
methods in a :ref:`execution-analysis-plan` as follows:

.. code-block:: python

    from mangadap.par.analysisplan import AnalysisPlanSet
    from mangadap.proc.spatiallybinnedspectra import SpatiallyBinnedSpectra
    from mangadap.proc.stellarcontinuummodel import StellarContinuumModel
    from mangadap.proc.emissionlinemodel import EmissionLineModel
    from mangadap.config.defaults import dap_method

    daptypes = []
    for plan in AnalysisPlanSet.from_par_file(plan_file):
        bin_method = SpatiallyBinnedSpectra.define_method(plan['bin_key'])
        sc_method = StellarContinuumModel.define_method(plan['continuum_key'])
        el_method = EmissionLineModel.define_method(plan['elfit_key'])
        daptypes += [dap_method(bin_method['key'], sc_method['template_library'],
                                el_method['continuum_templates'])]
