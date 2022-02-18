import os
import time
import numpy

import argparse

from scipy import interpolate

from mangadap.config import defaults
from mangadap.proc.templatelibrary import TemplateLibrary
from mangadap.util.sampling import spectral_coordinate_step

from mangadap.scripts import scriptbase


class TemplateFluxNorm(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Compute flux normalizations for DAP template set',
                                    width=width)
        parser.add_argument('key', type=str, help='DAP template-library key')
        parser.add_argument('ofile', type=str, help='Output file name')
        parser.add_argument('-v', '--velscale_ratio', default=4, type=int,
                            help='velocity scale ratio between the template pixel and '
                                 'the nominal pixel width')
        parser.add_argument('-s', '--step', default=1e-4, type=float,
                            help='linear or log-linear spectral pixel width')
        parser.add_argument('-l', '--linear', action='store_true', default=False,
                            help='pixel sampling should be linear (default is log)')
        return parser

    @staticmethod
    def main(args):

        # Get the g and r filter curves
        filter_root = os.path.join(defaults.dap_data_root(), 'filter_response')

        f = os.path.join(filter_root, 'gunn_2001_g_response.db')
        if not os.path.isfile(f):
            raise FileNotFoundError('No file: {0}'.format(f))
        db = numpy.genfromtxt(os.path.join(f))
        g = interpolate.interp1d(db[:,0], db[:,1], fill_value=0.0, bounds_error=False)

        f = os.path.join(filter_root, 'gunn_2001_r_response.db')
        if not os.path.isfile(f):
            raise FileNotFoundError('No file: {0}'.format(f))
        db = numpy.genfromtxt(os.path.join(f))
        r = interpolate.interp1d(db[:,0], db[:,1], fill_value=0.0, bounds_error=False)

        # Build the template library
        tpl = TemplateLibrary(args.key, velscale_ratio=args.velscale_ratio,
                              spectral_step=args.step, log=not args.linear, hardcopy=False)

        # Get the sampling
        nwave = len(tpl['WAVE'].data)
        dwave = spectral_coordinate_step(tpl['WAVE'].data, log=not args.linear)
        if not args.linear:
            dwave *= tpl['WAVE'].data*numpy.log(10.)

        flxmid = (numpy.amax(tpl['FLUX'].data, axis=1)+numpy.amin(tpl['FLUX'].data, axis=1))/2.0
        flxrng = numpy.amax(tpl['FLUX'].data, axis=1)-numpy.amin(tpl['FLUX'].data, axis=1)

        g_tc = g(tpl['WAVE'].data)
        r_tc = r(tpl['WAVE'].data)

        g_tc_sum = numpy.sum(g_tc*dwave)
        r_tc_sum = numpy.sum(r_tc*dwave)

        gflux = numpy.sum(tpl['FLUX'].data*g_tc[None,:]*dwave[None,:], axis=1)/g_tc_sum
        rflux = numpy.sum(tpl['FLUX'].data*r_tc[None,:]*dwave[None,:], axis=1)/r_tc_sum

        numpy.savetxt(args.ofile,
                      numpy.array([ numpy.arange(tpl.ntpl), flxmid, flxrng, gflux, rflux,
                                    -2.5*numpy.log10(gflux/rflux)]).T,
                      fmt=[ '%5d', '%15.8e', '%15.8e', '%15.8e', '%15.8e', '%15.8e'],
                      header='{0:>3} {1:>15} {2:>15} {3:>15} {4:>15} {5:>15}'.format(
                                'ID', 'FLXMID', 'FLXRNG', 'MEANg', 'MEANr', '-2.5LOG(g/r)'))

