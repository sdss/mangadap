import os
import time
import argparse

import numpy

from astropy.io import fits
import astropy.constants

from mangadap.datacube import MaNGADataCube
from mangadap.survey.drpcomplete import DRPComplete

from mangadap.scripts import scriptbase


class WriteDapConfig(scriptbase.ScriptBase):

    @classmethod
    def name(cls):
        """
        Return the name of the executable.
        """
        return 'write_dap_config'

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Generate a DAP input configuration file',
                                    width=width)

        parser.add_argument('plate', type=int, help='Plate number')
        parser.add_argument('ifudesign', type=int, help='IFU design number')
        parser.add_argument('ofile', type=str, help='Output file name')
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('-c', '--drpcomplete', type=str, help='DRP complete fits file',
                            default=None)
        group.add_argument('-a', '--drpall', type=str, help='DRPall fits file', default=None)

        parser.add_argument('--sres_ext', type=str, default=None,
                            help='Spectral resolution extension to use.  Default set by '
                                 'MaNGADataCube class.')
        parser.add_argument('--sres_fill', type=str, default=None,
                            help='If present, use interpolation to fill any masked pixels in the '
                                 'spectral resolution vectors. Default set by MaNGADataCube class.')
        parser.add_argument('--covar_ext', type=str, default=None,
                            help='Use this extension to define the spatial correlation matrix.  '
                                 'Default set by MaNGADataCube class.')
        parser.add_argument('--drpver', type=str, default=None,
                            help='DRP version.  Default set by MaNGADataCube class.')
        parser.add_argument('--redux_path', type=str, default=None,
                            help='Path to the top-level DRP reduction directory.  Default set by '
                                'MaNGADataCube class.')
        parser.add_argument('--directory_path', type=str, default=None,
                            help='Exact path to the directory with the MaNGA DRP datacube.  The '
                                 'name of the file itself must match the nominal MaNGA DRP naming '
                                 'convention.  Default set by MaNGADataCube class.')
        parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                            help='Overwrite any existing files.')
        return parser

    @staticmethod
    def main(args):
        t = time.perf_counter()

        if args.drpcomplete is not None:
            # Use the DRPcomplete file
            root_dir = os.path.dirname(args.drpcomplete)
            if len(root_dir) == 0:
                root_dir = '.'
            drpver = args.drpcomplete[args.drpcomplete.find('_v')+1:args.drpcomplete.find('.fits')]
            drpc = DRPComplete(drpver=drpver, directory_path=root_dir, readonly=True)
            index = drpc.entry_index(args.plate, args.ifudesign)
            MaNGADataCube.write_config(args.ofile, drpc['PLATE'][index], drpc['IFUDESIGN'][index],
                                       log=True,
                                       z=drpc['VEL'][index]/astropy.constants.c.to('km/s').value,
                                       vdisp=drpc['VDISP'][index], ell=drpc['ELL'][index],
                                       pa=drpc['PA'][index], reff=drpc['REFF'][index],
                                       sres_ext=args.sres_ext, sres_fill=args.sres_fill,
                                       covar_ext=args.covar_ext, drpver=args.drpver,
                                       redux_path=args.redux_path, overwrite=args.overwrite)
            return

        # Use the DRPall file
        with fits.open(args.drpall) as hdu:
            indx = numpy.where(hdu['MANGA'].data['PLATEIFU'] == '{0}-{1}'.format(args.plate,
                                                                                args.ifudesign))[0]
            if len(indx) != 1:
                raise ValueError('{0}-{1} either does not exist or has more than one match!'.format(
                                    args.plate, args.ifudesign))

            MaNGADataCube.write_config(args.ofile, args.plate, args.ifudesign,
                                       z=hdu[1].data['z'][indx[0]],
                                       ell=1-hdu[1].data['nsa_elpetro_ba'][indx[0]],
                                       pa=hdu[1].data['nsa_elpetro_phi'][indx[0]],
                                       reff=hdu[1].data['nsa_elpetro_th50_r'][indx[0]],
                                       sres_ext=args.sres_ext, sres_fill=args.sres_fill,
                                       covar_ext=args.covar_ext, drpver=args.drpver,
                                       redux_path=args.redux_path,
                                       directory_path=args.directory_path,
                                       overwrite=args.overwrite)

        print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))

