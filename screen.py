#!/usr/bin/env python3
import os
import sys
if not os.environ.get('MATSCREEN_PYTHON_TOKEN', False):
    SET_ENVIRONMENT = os.path.join(os.path.dirname(__file__), 'set_environment')
    exit(os.system("echo 'source %s; MATSCREEN_PYTHON_TOKEN=1 $MATSCREEN_PYTHON %s' | bash"%(SET_ENVIRONMENT, ' '.join(sys.argv))))

import argparse
import textwrap
from src.mpscreener import Screening

parser = argparse.ArgumentParser(
    description=textwrap.dedent('''\
         This is a python based command line tool to run screening simulations on a given set of processors available
         The are options to run grid, adsorption, or coadsorption simulations on a given set of structures with some chosen molecules. The number of steps of the simulation and the number of processes in a pool can be modified
         See example of usage in the run-python.sh file
         '''),
    epilog="")

parser.add_argument('-f', '--forcefield', nargs='?', default="UFF",
                    help='specify the force field you want to use (GenericMOFs, UFF, GenericZeolites, etc.)\nDefault = UFF')

parser.add_argument('-s', '--structures', nargs='?', required=True,
                    help='specify the structures to screen in a csv file\nExample: structures_sample.csv\nThis csv file is required and has to contain "Structures" as a header name.')

parser.add_argument('-n', '--nprocesses', nargs='?', required=True,
                    help='specify the number of processes to run at a time')

parser.add_argument('-ppn', '--procspernode', nargs='?', required=True,
                    help='specify the number of processors per node in the cluster')

parser.add_argument('-N', '--Ncycles', nargs='?', default='1000',
                    help='specify the number of cycles per GCMC calculation\nDefault = 1000')

parser.add_argument('-t', '--type', nargs='?', default='grid',
                    help='specify the type of simulation you want to run on the structures\nDefault = grid')

parser.add_argument('-m', '--molecules', nargs='+', default=['xenon','krypton'],
                    help='specify the adsorbent molecules\nExample: xenon krypton CO2.\nDefault=xenon krypton')

parser.add_argument('-T', '--Temperature', nargs='?', default=298.0,
                    help='specify the temperature for the Raspa2 simulations \nDefault=298.0')

parser.add_argument('-p', '--pressures', nargs='+', default=['101300'],
                    help='specify the pressures in the simulation\nExample: xenon krypton CO2. Default=101300')

parser.add_argument('-C', '--Cutoff', nargs='?', default=12,
                    help='specify the cut-off of the Raspa2 simulations\nDefault=12')

parser.add_argument('-c', '--composition', nargs=2, default=None,
                    help='specify the composition of each adsorbent molecule for coadsorption simulation\nExample: 90 10. Default=None')

parser.add_argument('-r', '--radius', nargs='?', default='1.2',
                    help='specify the radius of the probe in Zeo++ calculations or the sampling sphere relative radius\nDefault: 1.2')

parser.add_argument('-rj', '--rejection', nargs='?', default='0.85',
                    help='specify the rejection condition relative radius in surface simulations \nDefault: 0.85')

parser.add_argument('-R', '--restart', nargs='?', default='no',
                    help='specify if you want to restart from binary files')

parser.add_argument('-g', '--glost_list', nargs='?', default='no',
                    help='specify if you want to generate a glost list for mprun (cluster only)')

parser.add_argument('-o', '--output_directory', nargs='?', default='.',
                    help='specify the directory in which you want the simulation files to be installed. Default=. (current directory)')

args = parser.parse_args()

FORCE_FIELD = args.forcefield
structures_file = args.structures
option = args.type

MOLECULES = args.molecules
OUTPUT_PATH = args.output_directory
PRESSURES = args.pressures
temperature = float(args.Temperature)
cutoff = float(args.Cutoff)

CYCLES = int(args.Ncycles)
COMPOSITION = args.composition
radius = float(args.radius)
RESTART = args.restart

nprocs = int(args.nprocesses)
ppn = int(args.procspernode)
# initialising the screening procedure
screen = Screening(structures_file, ppn, nprocs, pressures=PRESSURES, temperature=temperature, cutoff=cutoff,probe_radius=radius,
         force_field=FORCE_FIELD, MOLECULES=MOLECULES, type_=option, composition=COMPOSITION, cycles=CYCLES, OUTPUT_PATH=OUTPUT_PATH,
         RESTART=RESTART)
# launching the screening
if args.glost_list in ['yes','y']:
  screen.glost_list()
else :
  screen.mp_run()

