import os
import time
import argparse

from mangadap.survey.drpcomplete import DRPComplete

def parse_args(options=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('drpcomplete', type=str, help='DRP complete fits file')
    parser.add_argument('plate', type=int, help='Plate number')
    parser.add_argument('ifudesign', type=int, help='IFU design number')
    parser.add_argument('mode', type=str, help='DRP 3D file mode: CUBE or RSS')
    parser.add_argument('ofile', type=str, help='Output file name')
    return parser.parse_args() if options is None else parser.parse_args(options)

def main(args):
    t = time.perf_counter()
    root_dir = os.path.dirname(args.drpcomplete)
    if len(root_dir) == 0:
        root_dir = '.'
    drpver = args.drpcomplete[args.drpcomplete.find('_v')+1 : args.drpcomplete.find('.fits')]
    drpc = DRPComplete(drpver=drpver, directory_path=root_dir, readonly=True)
    drpc.write_par(ofile=args.ofile, mode=args.mode, plate=args.plate, ifudesign=args.ifudesign)
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))

