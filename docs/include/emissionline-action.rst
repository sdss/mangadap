
The ``action`` parameter allows the emission-line database to be used
both in masking during the stellar-continuum modeling (see
:class:`~mangadap.util.pixelmask.SpectralPixelMask`) and during the
emission-line modeling itself.

The valid actions are:

    * ``i``: ignore the line, as if the line were commented out.
    * ``f``: fit the line and mask the line when fitting the stellar
      continuum.
    * ``m``: mask the line when fitting the stellar continuum but do
      *not* fit the line itself
    * ``s``: defines a sky line that should be masked.  When masked, the
      wavelength of the line is *not* adjusted for the redshift of the
      object spectrum.

I.e., when using the emission-line database for the emission-line
modeling, lines with the action set to ``f`` are fit, whereas all other
lines are ignored.

