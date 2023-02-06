
import warnings

from IPython import embed

import numpy
from scipy import sparse, linalg

from mangadap.scripts import scriptbase


def impose_positive_definite(mat, gpm=None, min_eigenvalue=1e-10, renormalize=True, maxiter=1,
                             quiet=False):
    """
    Force a matrix to be positive-semidefinite.

    Following, e.g.,
    http://comisef.wikidot.com/tutorial:repairingcorrelation, the algorithm
    is as follows:

        - Calculate the eigenvalues and eigenvectors of the provided matrix
          (this is the most expensive step).
        - Impose a minimum eigenvalue (see ``min_eigenvalue``)
        - Reconstruct the input matrix using the eigenvectors and the
          adjusted eigenvalues
        - Renormalize the reconstructed matrix such that its diagonal is
          identical to the input matrix, if requested.
        - Iterate on the above until the adjusted matrix is
          positive-semidefinite or for the provided maximum number of
          iterations.

    Args:
        mat (`numpy.ndarray`_, `scipy.sparse.csr_matrix`_):
            The matrix to force to be positive definite.
        gpm (`numpy.ndarray`_, optional):
            A boolean array selecting the entries along the diagonal to
            consider.  The code ignores the rows and columns associated
            with the diagonal elements that are *not* selected by this
            vector.  The adjustments to make the matrix
            positive-semidefinite is then only applied to the relevant
            submatrix.  This is useful for selecting non-zero elements
            of the diagonal.
        min_eigenvalue (:obj:`float`, optional):
            The minimum allowed matrix eigenvalue.
        renormalize (:obj:`bool`, optional):
            Include the renormalization (last) step in the list above.
        maxiter (:obj:`int`, optional):
            The maximum number of iterations to perform.
        quiet (:obj:`bool`, optional):
            Suppress output messages

    Returns:
        `numpy.ndarray`_, `scipy.sparse.csr_matrix`_: The modified
        matrix.  The return time should match the input type.
    """
    _mat = mat if gpm is None else mat[numpy.ix_(gpm,gpm)]
    _inp_diag = _mat.diagonal()
    # TODO: Check the input types?
    if is_positive_definite(_mat):
        return mat
    _mat = _mat.toarray()
    i = 0
    while i < maxiter and not is_positive_definite(_mat):
        # Get the eigenvalues/eigenvectors.  NOTE: This command can take
        # a while, depending on the size of the array...
        w, v = numpy.linalg.eig(_mat)
        if not quiet:
            # Report the number of non-positive values
            indx = numpy.logical_not(numpy.real(w) > 0)
            print(f'Iteration {i+1} found {numpy.sum(indx)} non-positive eigenvalues.')
        if numpy.all(numpy.real(w) > 0):
            break
        # Clip to a minimum eigenvalue
        w = numpy.maximum(w, min_eigenvalue)
        # Reconstruct with the new eigenvalues, keeping only the real
        # component...
        _mat = numpy.real(numpy.dot(v, numpy.dot(numpy.diag(w), v.T)))
        if renormalize:
            # Renormalize the matrix so that the diagonals are identical
            r = numpy.sqrt(_inp_diag/_mat.diagonal())
            _mat *= numpy.outer(r,r)
        i += 1

    if gpm is None:
        return sparse.csr_matrix(_mat) if sparse.issparse(mat) else _mat

    pdmat = mat.copy()
    pdmat[numpy.ix_(gpm,gpm)] = _mat
    return pdmat


def is_positive_definite(mat, quiet=True, quick=True):
    r"""
    Check if a matrix is positive definite.

    This is done by calculating the eigenvalues and eigenvectors of the
    provided matrix and checking if all the eigenvalues are :math:`>0`.
    Because of that, it is nearly as expensive as just calling
    :func:`impose_positive_definite`.

    Args:
        mat (`numpy.ndarray`_, `scipy.sparse.csr_matrix`_):
            The matrix to check.
        quiet (:obj:`bool`, optional):
            Suppress terminal output.
        quick (:obj:`bool`, optional):
            Use the quick method, which is to try to use Cholesky decomposition
            and check if it throws a LinAlgError.  The slow way is to determine
            the eigenvalues and check if they are all positive.  If True and
            quiet is False, only the error reported by the Cholesky
            decomposition is printed, instead of the full list of non-positive
            eigenvalues.

    Returns:
        :obj:`bool`: Flag that matrix is positive definite.
    """
    _mat = mat.toarray() if isinstance(mat, sparse.csr_matrix) else mat

    if quick:
        try:
            cho = linalg.cholesky(_mat)
        except linalg.LinAlgError as e:
            if not quiet:
                print(str(e))
            return False
        else:
            return True

    # Get the eigenvalues/eigenvectors
    w, v = numpy.linalg.eig(_mat)
    notpos = numpy.logical_not(numpy.real(w) > 0)
    if not quiet:
        if numpy.any(notpos):
            warnings.warn(f'{numpy.sum(notpos)} eigenvalues are not positive!')
            print('{0:>6} {1:>8}'.format('Index', 'EigenVal'))
            for i in numpy.where(notpos)[0]:
                print('{0:>6} {1:8.2e}'.format(i, w[i]))
    return not numpy.any(notpos)


class MangaSynthDatacube(scriptbase.ScriptBase):

    @classmethod
    def name(cls):
        """
        Return the name of the executable.
        """
        return 'manga_synth_datacube'

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Create a synthetic datacube', width=width)

        parser.add_argument('plateifu', help='PLATE-IFU number')
        parser.add_argument('oroot', help='Root name for output files')
        parser.add_argument('-d', '--directory_path', type=str, default=None,
                            help='Path to CUBE and RSS file')
        parser.add_argument('-e', '--error', help='Multiplicative factor to apply to error',
                            default=1., type=float)
        parser.add_argument('-n', '--nsim', help='Number of realizations to make', default=1,
                            type=int)

        return parser

    @staticmethod
    def main(args):

        from pathlib import Path

        import numpy
        from scipy import sparse, spatial, interpolate

        from astropy.io import fits

        from mangadap.datacube import MaNGADataCube
        from mangadap.spectra import MaNGARSS
        from mangadap.proc.templatelibrary import TemplateLibrary
        from mangadap.proc.ppxffit import PPXFFit
        from mangadap.proc.spectralfitting import StellarKinematicsFit
        from mangadap.util.geometry import projected_polar
        from mangadap.util.sampling import spectral_coordinate_step
        from mangadap.util.fileio import compress_file

        odir = Path(args.oroot).resolve()
        oroot = odir.name
        odir = odir.parent
        if not odir.exists():
            odir.mkdir(parents=True)

        # Maximum size for the random draw array
        max_gib = 50.
#        max_gib = 0.1

        # Read the cube 
        plate, ifu = map(lambda x : int(x), args.plateifu.split('-'))
        cube = MaNGADataCube.from_plateifu(plate, ifu, directory_path=args.directory_path)

        # Get the cube on-sky coordinates (relative to the object center)
        cube_x, cube_y = cube.mean_sky_coordinates()
        # And a good-pixel mask for the spaxels with *any* valid data
#        spat_gpm = numpy.any(numpy.logical_not(cube.mask > 0), axis=-1)
        spat_gpm = numpy.any(cube.flux != 0, axis=-1)

        # Read the template library
        velscale_ratio = 2
        nwave = cube.wave.size
        obj_wave = cube.wave
        tpl = TemplateLibrary('MASTARHC2', velscale_ratio=velscale_ratio,
                              spectral_step=spectral_coordinate_step(obj_wave, log=cube.log),
                              log=True, hardcopy=False)

        # Build the synthetic cube based on a single spectrum
        tpl_indx = 39
        tpl_wave = tpl['WAVE'].data
        tpl_flux = tpl['FLUX'].data[tpl_indx].reshape(1,-1)

        # Set the true velocity and velocity-dispersion field
        inc = 40. # deg
        pa = 45. # deg
        vmax = 100. # km/s
        disp = 100. # km/s

        # Use the cube coordinates to build the expected kinematics
        spat_ij = numpy.ravel_multi_index(numpy.where(spat_gpm), spat_gpm.shape)
        xref = cube_x[spat_gpm]
        yref = cube_y[spat_gpm]
        ngpm = numpy.sum(spat_gpm)

        r, th = projected_polar(xref, yref, *numpy.radians([pa, inc]))
        maxr = numpy.amax(r)
#        v = r*vmax*numpy.cos(th)/maxr
        v = numpy.full(r.shape, 300.)
        disp = numpy.full(v.shape, disp, dtype=float)
        v_map = numpy.zeros(spat_gpm.shape, dtype=float)
        v_map[spat_gpm] = v
        disp_map = numpy.zeros(spat_gpm.shape, dtype=float)
        disp_map[spat_gpm] = disp

        model_par = StellarKinematicsFit.init_datatable(1, 0, 0, 2, numpy.int16, shape=ngpm)
        model_par['BINID'] = spat_ij

        cube_flux = numpy.ma.MaskedArray(cube.flux.reshape(-1,4563))[spat_gpm.ravel(),:]
        cube_flux[numpy.logical_not(cube_flux != 0)] = numpy.ma.masked
#        model_par['TPLWGT'] = numpy.ma.median(cube_flux, axis=1).filled(0.0).reshape(-1,1)
        model_par['TPLWGT'] = 1.
        flux_map = numpy.zeros(spat_gpm.shape, dtype=float)
        flux_map[spat_gpm] = model_par['TPLWGT'][:,0]

        gpm = numpy.logical_not(numpy.ma.getmaskarray(cube_flux))
        ntpl_pix = (tpl_wave.size - tpl_wave.size % velscale_ratio) // velscale_ratio
        for i in range(ngpm):
            model_par['BEGPIX'][i] = numpy.where(gpm[i])[0][0]
            model_par['ENDPIX'][i] = numpy.where(gpm[i])[0][-1]
            while ntpl_pix < model_par['ENDPIX'][i] - model_par['BEGPIX'][i]:
                model_par['BEGPIX'][i] += 15
                model_par['ENDPIX'][i] -= 15
        model_par['BEGPIX'] += 100
        model_par['ENDPIX'] -= 100
        model_par['KIN'] = numpy.column_stack((v, disp))

        cube_models = PPXFFit.construct_models(tpl_wave, tpl_flux, obj_wave, cube_flux.shape,
                                    model_par, select=model_par['BEGPIX'] < model_par['ENDPIX'])
        cube_flux = numpy.ma.masked_all(cube.flux.shape, dtype=float).reshape(-1, 4563)
        cube_flux[spat_ij] = cube_models
        cube_flux = cube_flux.reshape(cube.flux.shape)

        # Read the RSS files
        cube.load_rss()

        # Get the error in the rss spectra needed to create the covariance
        # matrix in the datacubes
        obj_err = args.error * numpy.sqrt(numpy.ma.divide(1, cube.rss.ivar).filled(0.0))

        # Instantiate the random number generator
        rng = numpy.random.default_rng()

        # Instantiate the output arrays
        spatial_shape = cube.spatial_shape
        cube_shape = spatial_shape + (nwave,)

        # One cube per simulations
#        flux = numpy.zeros((args.nsim,)+cube_shape, dtype=numpy.float32)
        flux = cube_flux.filled(0.0).astype(numpy.float32).reshape(-1,4563)

        # The variance and mask arrays are identical for all simulations
#        var = numpy.zeros(cube_shape, dtype=numpy.float32).reshape(-1,4563)
        var = numpy.zeros(cube_shape, dtype=numpy.float32)
        mask = numpy.ma.getmaskarray(cube_flux).copy()
#        bad_draw = numpy.zeros(nwave, dtype=bool)

        nsim = numpy.array([args.nsim])
        # Size of a float64 in GiB
        float64_size = numpy.dtype(numpy.float64).itemsize/2**30
        while numpy.prod((nsim[0],) + cube.rss.shape) * float64_size > max_gib:
            nsim = numpy.array([[n//2,n//2 + n%2] for n in nsim]).ravel()

        n_written = 0
        for k in range(nsim.size):
            print(f'Working on subset {k+1} of {nsim.size}.')

            _flux = numpy.tile(flux, (nsim[k],1,1))
            draw = rng.normal(size=(nsim[k],) + cube.rss.shape)
            draw *= obj_err[None,...]

            # Iterate over wavelength channels
            for j in range(nwave):
                print(f'Wave: {j+1}/{nwave}', end='\r')
                if numpy.all(mask[...,j]):
                    continue
                # Get the rectification matrix
                t = cube.rss.rectification_transfer_matrix(j, quiet=True)
                # Get the covariance in the cube for this channel
                covar = t.dot(sparse.diags(obj_err[:,j]**2).dot(t.T))

                # ------------------------------------------------------------------
                # NEW APPROACH
                var[...,j] = covar.diagonal().reshape(spatial_shape)
                for i in range(nsim[k]):
                    _flux[i,:,j] += t.dot(draw[i,:,j])
                # ------------------------------------------------------------------

            print(f'Wave: {nwave}/{nwave}')

            # Reshape the flux array
            _flux = _flux.reshape((nsim[k],) + cube_shape)
            _flux[:,mask] = 0.
    #        var = var.reshape(cube_shape)

            # Copy the WCS, wave, sres
            ivar = numpy.ma.divide(1, var).filled(0.0)
            for i in range(nsim[k]):
                ofile = odir / f'{oroot}-{n_written+1:02}.fits'
                fits.HDUList([fits.PrimaryHDU(),
                            fits.ImageHDU(data=obj_wave, name='WAVE'),
                            fits.ImageHDU(data=_flux[i], name='FLUX', header=cube.wcs.to_header()),
                            fits.ImageHDU(data=ivar, name='IVAR'),
                            fits.ImageHDU(data=mask.astype(numpy.int16), name='MASK'),
                            fits.ImageHDU(data=cube.sres.astype(numpy.float32), name='SRES'),
    #                          fits.ImageHDU(data=bad_draw.astype(numpy.int16), name='DRAW'),
                            fits.ImageHDU(data=flux_map.astype(numpy.float32), name='INP_WGT'),
                            fits.ImageHDU(data=v_map.astype(numpy.float32), name='INP_V'),
                            fits.ImageHDU(data=disp_map.astype(numpy.float32), name='INP_SIG')
                            ]).writeto(str(ofile), overwrite=True, checksum=True)
                compress_file(str(ofile), overwrite=True, rm_original=True)            
                n_written += 1


