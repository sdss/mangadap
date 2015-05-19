#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
This script is used to launch dapqa.

"""

import jinja2
import os

templateLoader = jinja2.FileSystemLoader(os.path.join(os.getcwd(), 'templates'))
templateEnv = jinja2.Environment( loader=templateLoader )

TEMPLATE_FILE = 'maps.jinja2'
template = templateEnv.get_template( TEMPLATE_FILE )

# Here we add a new input variable containing a list.
# Its contents will be expanded in the HTML as a unordered list.

plates = [7443]
ifus = [6102, 12702]
dap_modes = ['CUBE']
binnings = ['NONE-002', 'STON-001']

plates = map(str, plates)
ifus = map(str, ifus)

topdir = os.getenv('MANGA_SPECTRO_ANALYSIS')

plot_types = [
    'snr_maps.png',
    'kin_maps.png',
    'emflux_ew_maps.png',
    'emflux_fb_maps.png']




for plate in plates:
    for ifu in ifus:
        path_plots = os.path.join(topdir, plate, ifu, 'plots', '')
        for dap_mode in dap_modes:
            for binning in binnings:
                stem = ('-'.join(['manga', plate, ifu,
                        'LOG%s_BIN' % dap_mode, binning]) + '_')
                urls = []
                for plot_type in plot_types:
                    urls.append(path_plots + stem + plot_type)
                templateVars = dict(urls=urls)
                outputText = template.render(templateVars)
                f = open(path_plots + binning + '_maps.html', 'w')
                f.write(outputText)
                f.close()


# templateVars = dict(urls=urls)
# 
# 
# outputText = template.render( templateVars )
# 
# f = open(path_plots + 'maps.html', 'w')
# f.write(outputText)
# f.close()
