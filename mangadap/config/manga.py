"""
Class defining the I/O configuration for MaNGA files in the DAP.
"""
from pathlib import Path
from astropy.io import fits

from . import defaults
from ..util.parser import DefaultConfig

class MaNGAConfig:

    instrument = 'manga'

    def __init__(self, plate, ifudesign, log=True, drpver=None, redux_path=None,
                 directory_path=None, analysis_path=None):
        """
        Initialize the I/O configuration using the same file as used to read
        MaNGA data.
        """
        # MaNGA-specific attributes
        self.plate = plate
        self.ifudesign = ifudesign
        self.log = log
        self.drpver = defaults.drp_version() if drpver is None else drpver
        self.directory_path = defaults.drp_directory_path(self.plate, drpver=self.drpver,
                                                          redux_path=redux_path) \
                                if directory_path is None else Path(directory_path).resolve()
        self.analysis_path = defaults.dap_analysis_path(drpver=self.drpver) \
                                if analysis_path is None else Path(analysis_path).resolve()
        self.mode = 'CUBE'
        self.ebv_key = 'EBVGAL'

    @classmethod
    def from_config(cls, cfgfile):
        """
        """
        # Read the configuration file
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
        return cls(plate, ifu, log=log, drpver=drpver, redux_path=redux_path,
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

        if 'CUBE' not in datafile:
            raise NotImplementedError('Currently can only analyze MaNGA CUBE files.')

        # Parse the relevant information from the filename
        plate, ifudesign, log = ifile.name.split('-')[1:4]

        if 'LOG' not in log:
            raise NotImplementedError('Currently can only analyze LOG binned files.')

        with fits.open(str(ifile)) as hdu:
            drpver = hdu[0].header['VERSDRP3']
        # Instantiate the DRPFits base
        return cls(int(plate), int(ifudesign), log='LOG' in log, drpver=drpver,
                   directory_path=ifile.parent)

    @property
    def file_name(self):
        MaNGAConfig.check_mode(self.mode)
        root = defaults.manga_fits_root(self.plate, self.ifudesign,
                                        mode=f'{"LOG" if self.log else "LIN"}{self.mode}')
        return f'{root}.fits.gz'

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
    def identifier(self):
        return f'{self.instrument}-{self.key}'


