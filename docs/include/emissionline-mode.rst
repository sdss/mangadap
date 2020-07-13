
The ``mode`` parameter sets how the emission line should be treated with
respect of the rest of the lines being modeled.

The valid modes are:

    * ``f``: Fit the line independently of all others.
    * ``wN``: Used by :class:`~mangadap.proc.elric.Elric` only.  Fit the
      line with untied parameters, but use a window that includes both
      this line and the line with index ``N``.
    * ``xN``: Used by :class:`~mangadap.proc.elric.Elric` only.  Fit the
      line with its flux tied to the line with index ``N``.
    * ``vN``: Fit the line with the velocity tied to the line with index
      ``N``.
    * ``sN``: Fit the line with the velocity dispersion tied to the line
      with index ``N``.
    * ``kN``: Fit the line with the velocity and velocity dispersion
      tied to the line with index ``N``.
    * ``aN``: Fit the line with the flux, velocity, and velocity
      dispersion tied to the line with index ``N``.

As noted in the mode description, many of the modes are only
available when using the :class:`~mangadap.proc.elric.Elric` module.
For the ``w`` mode, this is simply because the preferred module,
:class:`~mangadap.proc.sasuke.Sasuke`, fits the full spectrum instead
of fitting the lines within small spectral windows. The other
limitation are because :class:`~mangadap.proc.sasuke.Sasuke` is based
on the use of template spectra to fit the emission lines (see
:class:`~mangadap.proc.emissionlinetemplates.EmissionLineTemplates`):
To tie line fluxes, the lines to be tied are included in the same
template spectrum, meaning that their kinematics are also
automatically tied. That means that, for
:class:`~mangadap.proc.sasuke.Sasuke`, the ``x`` and ``a`` modes are
identical.

In the :ref:`emissionlines-model-par-format` example, the modes set
the :math:`{\rm H}\alpha` line as the "reference" line. I.e., there
should always be one line whose mode is ``f``. This requirement is
simply practical in setting up the tied parameter structure; there is
no more weight given to the fit to the reference line than any other
line. The blue [OII] and red [NII] lines have their velocities tied
to the :math:`{\rm H}\alpha` line, all kinematics of the red [OII]
line are tied to the blue [OII] line, and all parameters of the blue
[NII] line are tied to the red [NII] line with a fixed flux ratio of
``NII-6550/NII-6585 == 0.34``. By virtue of being tied to lines that
have their velocites tied to the :math:`{\rm H}\alpha` line, the
velocities of the red [OII] and blue [NII] lines are also tied to the
:math:`{\rm H}\alpha` line. Again, this doesn't mean that the fit to
the :math:`{\rm H}\alpha` line is given any more weight than any
other line, it just means that there is one model parameter that
defines the velocity of *all* lines.

