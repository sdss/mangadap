
from ..util.pkg import all_subclasses
from . import scriptbase

# The import of all the script modules here is what enables the dynamic
# compiling of all the available scripts below
from . import calculate_covariance
from . import construct_dapall
from . import dap_status
from . import dapall_qa
from . import find_repeat_observations
from . import fit_residuals
from . import manga_dap_inspector
from . import manga_dap
from . import plate_fit_qa
from . import ppxffit_qa
#from . import rundap
from . import spotcheck_dap_maps
from . import template_flux_norm
from . import write_dap_config

# Build the list of script classes
def script_classes():
    import numpy as np

    # Recursively collect all subclasses
    scr_c = np.array(list(all_subclasses(scriptbase.ScriptBase)))
    scr_n = np.array([c.name() for c in scr_c])
    # Construct a dictionary with the script name and class
    srt = np.argsort(scr_n)
    return dict([ (n,c) for n,c in zip(scr_n[srt],scr_c[srt])])

mangadap_scripts = list(script_classes().keys())

