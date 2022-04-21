
from IPython import embed

from mangadap.scripts import scriptbase


# TODO: Move this somewhere that can be more generally accessed?
def calculate_covariance_cube(plate, ifudesign, ofile, nc=1, wave=None, directory_path=None):

    import numpy
    from mangadap.spectra import MaNGARSS

    print('     PLATE: {0}'.format(plate))
    print(' IFUDESIGN: {0}'.format(ifudesign))

    # Access the DRP RSS file
    print('Attempting to open RSS file:')
    rss = MaNGARSS.from_plateifu(plate, ifudesign, directory_path=directory_path)
    print('     FOUND: {0}'.format(rss.file_path()))

    if wave is not None:
        channel = numpy.argsort(numpy.absolute(rss.wave - wave))[0]
        print('Nearest wavelength channel has wavelength {0:.1f} ang.'.format(rss.wave[channel]))
        C = rss.covariance_matrix(channel)
    else:
        if nc >= rss.nwave or nc == 0:
            print('Calculating full covariance cube ...')
            C = rss.covariance_cube()
            print('... done.')
        elif nc == 1:
            channel = rss.nwave//2
            print('Calculating covariance matrix at central channel: '
                  '{0:.1f} ang.'.format(rss.wave[channel]))
            C = rss.covariance_matrix(channel)
        else:
            print('Calculating covariance in {0} wavelength channels...'.format(nc))
            channels = numpy.linspace(0, rss.nwave-1, num=nc, dtype=int)
            C = rss.covariance_cube(channels=channels)

    print('Writing data to {0}.'.format(ofile))
    C.write(ofile, clobber=True)            # Write the data


class CalculateCovariance(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(description='Calculate datacube covariance', width=width)

        parser.add_argument('plate', type=int, help='plate ID to process')
        parser.add_argument('ifudesign', type=int, help='IFU design to process')
        parser.add_argument('output_file', type=str, help='Name for output file')

        mode = parser.add_mutually_exclusive_group(required=False)
        mode.add_argument('-n', '--numchannels', type=int,
                        help='Number of channels spread across the wavelength range for which '
                            'to compute the covariance matrix.  A value of 0 forces construction '
                            'of the full covariance cube.  The default is to calculate the '
                            'covariance matrix for a single channel at the central wavelength',
                        default=1)
        mode.add_argument('-w', '--wavelength', type=float,
                        help='Wavelength at which to compute a single covariance matrix; default '
                            'is the central wavelength', default=None)

        parser.add_argument('-d', '--directory_path', type=str, default=None,
                            help='Directory with the DRP produced RSS file; default uses '
                                 'environmental variables to define the default MaNGA DRP '
                                 'redux path')
        return parser

    @staticmethod
    def main(args):

        import os
        import time

        t = time.perf_counter()
        if os.path.isfile(args.output_file):
            print('WARNING: Overwriting existing file {0}!'.format(args.output_file))
        calculate_covariance_cube(args.plate, args.ifudesign, args.output_file, nc=args.numchannels,
                                wave=args.wavelength, directory_path=args.directory_path)
        print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))


