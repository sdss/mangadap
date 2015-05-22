"""
Create html pages that flip through maps of one galaxy.

"""
from __future__ import print_function

import jinja2
import os
import json

dict_tmp = {}
try:
    dict_tmp.iteritems
except AttributeError:
    def iterkeys(d):
        return iter(d.keys())
else:
    def iterkeys(d):
        return d.iterkeys()

templateLoader = jinja2.FileSystemLoader(os.path.join(os.getcwd(), 'templates'))
templateEnv = jinja2.Environment( loader=templateLoader )

TEMPLATE_FILE = 'maps.jinja2'
template = templateEnv.get_template( TEMPLATE_FILE )

path_data = os.path.join(os.getenv('MANGA_SPECTRO_ANALYSIS'), 'trunk_mpl3')
# topdir_sas = path_data
topdir_sas = 'https://data.sdss.org/sas/' + path_data.split('sdss03')[1]

plot_types = [
    'snr_maps',
    'kin_maps',
    'emflux_ew_maps',
    'emflux_fb_maps']
grad_types = ['emflux_gradients']
global_spec = ['spec-0000']

plot_type_filenames = [it + '.png' for it in plot_types]
grad_type_filenames = [it + '.png' for it in grad_types]
global_spec_filenames = [it + '.png' for it in global_spec]


#data = [{'plate':7443, 'ifudsgn':'12702', 'mode':'CUBE', 'binning':'NONE-002', 'spec':'0000', 'plotname':'maps'},
#        {'plate':7443, 'ifudsgn':'6102', 'mode':'CUBE', 'binning':'STON-001', 'spec':'0000', 'plotname':'maps'},
#        {'plate':7443, 'ifudsgn':'6102', 'mode':'CUBE', 'binning':'NONE-002', 'spec':'0000', 'plotname':'maps'}, 
#        {'plate':7443, 'ifudsgn':'3701', 'mode':'CUBE', 'binning':'RADIAL-003', 'spec':'0000', 'plotname':'maps'},
#        {'plate':7443, 'ifudsgn':'3701', 'mode':'CUBE', 'binning':'RADIAL-004', 'spec':'0000', 'plotname':'maps'}, 
#        {'plate':7443, 'ifudsgn':'3701', 'mode':'CUBE', 'binning':'ALL-005', 'spec':'0000', 'plotname':'maps'},
#        {'plate':7443, 'ifudsgn':'3701', 'mode':'CUBE', 'binning':'ALL-006', 'spec':'0000', 'plotname':'maps'},
#        {'plate':7443, 'ifudsgn':'3701', 'mode':'CUBE', 'binning':'ALL-007', 'spec':'0000', 'plotname':'maps'},
#        {'plate':7443, 'ifudsgn':'3701', 'mode':'RSS', 'binning':'RADIAL-001', 'spec':'0000', 'plotname':'maps'},
#        {'plate':7443, 'ifudsgn':'3701', 'mode':'RSS', 'binning':'RADIAL-002', 'spec':'0000', 'plotname':'maps'}, 
#        {'plate':7443, 'ifudsgn':'3701', 'mode':'RSS', 'binning':'ALL-003', 'spec':'0000', 'plotname':'maps'},
#        {'plate':7443, 'ifudsgn':'3701', 'mode':'RSS', 'binning':'ALL-004', 'spec':'0000', 'plotname':'maps'},
#        {'plate':7443, 'ifudsgn':'3701', 'mode':'RSS', 'binning':'ALL-005', 'spec':'0000', 'plotname':'maps'},
#        ]

fin = open(os.path.join(path_data, 'table_data.json'), 'r')
data = json.load(fin)
fin.close()




cube_ptypes = {
    'STON-001': plot_types,
    'NONE-002': plot_types, 
    'RADIAL-003': grad_types,
    'RADIAL-004': grad_types,
    'ALL-005': global_spec,
    'ALL-006': global_spec,
    'ALL-007': global_spec
    }

rss_ptypes = {
    'RADIAL-001': grad_types,
    'RADIAL-002': grad_types,
    'ALL-003': global_spec,
    'ALL-004': global_spec,
    'ALL-005': global_spec
    }



all_urls = {
    'CUBE': {
        'STON-001': {},
        'NONE-002': {}, 
        'RADIAL-003': {},
        'RADIAL-004': {},
        'ALL-005': {},
        'ALL-006': {},
        'ALL-007': {}
        },
    'RSS': {
        'RADIAL-001': {},
        'RADIAL-002': {},
        'ALL-003': {},
        'ALL-004': {},
        'ALL-005': {}
        }
    }

for p in plot_types:
    for binning in ['STON-001', 'NONE-002']:
        all_urls['CUBE'][binning][p] = []

for g in grad_types:
    for binning in ['RADIAL-003', 'RADIAL-004']:
        all_urls['CUBE'][binning][g] = []
    for binning in ['RADIAL-001', 'RADIAL-002']:
        all_urls['RSS'][binning][g] = []

for binning in ['ALL-005', 'ALL-006', 'ALL-007']:
    all_urls['CUBE'][binning]['spec-0000'] = []

for binning in ['ALL-003', 'ALL-004', 'ALL-005']:
    all_urls['RSS'][binning]['spec-0000'] = []



for it in data:
    path_plots_utah = os.path.join(path_data, str(it['plate']), it['ifudsgn'], 'plots', '')
    path_plots_sas = os.path.join(topdir_sas, str(it['plate']), it['ifudsgn'], 'plots', '')
    stem = ('-'.join(['manga', str(it['plate']), it['ifudsgn'],
                        'LOG%s_BIN' % it['mode'], it['binning']]) + '_')
    if it['mode'] == 'CUBE':
        ptypes = cube_ptypes[it['binning']]
    elif it['mode'] == 'RSS':
        ptypes = rss_ptypes[it['binning']]
    urls = []
    for i, plot_type in enumerate(ptypes):
        specdir = ''
        if 'ALL' in it['binning']:
            specdir = 'spectra/'
        url = path_plots_sas + specdir + stem + plot_type + '.png'
        urls.append(url)
        all_urls[it['mode']][it['binning']][plot_type].append(url)
    templateVars = dict(urls=urls)
    outputText = template.render(templateVars)
    if it['binning'] in ['STON-001', 'NONE-002']:
        suffix = '_maps'
    elif 'RADIAL' in it['binning']:
        suffix = '_gradients'
    else:
        suffix = '_global_spec'
    f = open(path_plots_utah + it['binning'] + suffix + '.html', 'w')
    f.write(outputText)
    f.close()


for k1 in iterkeys(all_urls):
    for k2 in iterkeys(all_urls[k1]):
        for k3 in iterkeys(all_urls[k1][k2]):
            templateVars = dict(urls=all_urls[k1][k2][k3])
            outputText = template.render(templateVars)
            f = open(os.path.join(path_data, k1+'_allgal_plots', '_'.join([k2, k3])) + '.html', 'w')
            f.write(outputText)
            f.close()

