"""
Class defining the I/O configuration for MaNGA files in the DAP.
"""
import warnings
from pathlib import Path
from astropy.io import fits

from . import defaults
from ..util.parser import DefaultConfig
from ..util.drpbitmask import DRPQuality3DBitMask

class MaNGAConfig:

    instrument = 'manga'

    def __init__(self, plate, ifudesign, mode='CUBE', log=True, drpver=None, redux_path=None,
                 chart_path=None, directory_path=None, analysis_path=None):
        """
        Initialize the I/O configuration using the same file as used to read
        MaNGA data.
        """
        # MaNGA-specific attributes
        self.plate = plate
        self.ifudesign = ifudesign
        MaNGAConfig.check_mode(mode)
        self.mode = mode
        self.log = log
        self.drpver = defaults.drp_version() if drpver is None else drpver
        # TODO: Depending on what is provided to the function, redux path, chart
        # path, and directory path can all be inconsistent.
        self.redux_path = defaults.drp_redux_path(drpver=self.drpver) \
                            if redux_path is None else Path(redux_path).resolve()
        self.chart_path = defaults.drp_finding_chart_path(self.plate, drpver=self.drpver,
                                                          redux_path=self.redux_path) \
                                if chart_path is None else Path(chart_path).resolve()
        self.directory_path = defaults.drp_directory_path(self.plate, drpver=self.drpver,
                                                          redux_path=self.redux_path) \
                                if directory_path is None else Path(directory_path).resolve()
        # TODO: This is only relevant the CUBE data
        if self.directory_path != defaults.drp_directory_path(self.plate, drpver=self.drpver,
                                                              redux_path=redux_path):
            warnings.warn('Paths do not match MaNGA defaults.  May not be able to find paired '
                          'CUBE & RSS files, if requested.')

        self.analysis_path = defaults.dap_analysis_path(drpver=self.drpver) \
                                if analysis_path is None else Path(analysis_path).resolve()
        # TODO: Turn these all into meta
        self.ebv_key = 'EBVGAL'
        self.qual_bitmask = DRPQuality3DBitMask()
        self.qual_key = 'DRP3QUAL'
        self.qual_flag = 'CRITICAL'


    @classmethod
    def from_config(cls, cfgfile, mode='CUBE'):
        """
        """
        # Read the configuration file.  This checks if the file exists and will
        # raise an exception if it does not.
        cfg = DefaultConfig(cfgfile, interpolate=True)

        # Set the attributes, forcing a known type
        plate = cfg.getint('plate')
        ifu = cfg.getint('ifu')
        if plate is None or ifu is None:
            raise ValueError('Configuration file must define the plate and IFU.')
        log = cfg.getbool('log', default=True)
        drpver = cfg.get('drpver')
        redux_path = cfg.get('redux_path')
        directory_path = cfg.get('directory_path')
        return cls(plate, ifu, mode=mode, log=log, drpver=drpver, redux_path=redux_path,
                   directory_path=directory_path)

    @classmethod
    def from_file(cls, datafile):
        """
        Initialize the I/O configuration from the datafile for the DAP to
        process.
        """
        ifile = Path(datafile).resolve()
        if not ifile.exists():
            raise FileNotFoundError(f'File does not exist: {ifile}')

        # Parse the relevant information from the filename
        plate, ifudesign, log = ifile.name.split('-')[1:4]

        with fits.open(str(ifile)) as hdu:
            drpver = hdu[0].header['VERSDRP3']

        # Instantiate the DRPFits base
        return cls(int(plate), int(ifudesign), mode='CUBE' if 'CUBE' in log else 'RSS',
                   log='LOG' in log, drpver=drpver, directory_path=ifile.parent)

    def copy(self):
        """
        Create a deep copy.
        """
        return MaNGAConfig(self.plate, self.ifudesign, log=self.log, drpver=self.drpver,
                           directory_path=self.directory_path, analysis_path=self.analysis_path)

    @property
    def file_name(self):
        root = defaults.manga_fits_root(self.plate, self.ifudesign,
                                        mode=f'{"LOG" if self.log else "LIN"}{self.mode}')
        return f'{root}.fits.gz'

    @property
    def file_path(self):
        return self.directory_path / self.file_name

    @staticmethod
    def check_mode(mode):
        """
        Check that the mode is valid.

        Valide modes are set by :func:`mode_options`.

        Args:
            mode (:obj:`str`):
                Mode value to check.

        Raises:
            ValueError:
                Raised if the mode is undefined.
        """
        options = ['CUBE', 'RSS']
        if mode not in options:
            raise ValueError(f'Unknown mode {mode}.  Must be in: {options}')

    @property
    def key(self):
        """
        Unique key used to identify the input data.
        """
        return f'{self.plate}-{self.ifudesign}'

    @property
    def output_root(self):
        return f'{self.instrument}-{self.key}'


