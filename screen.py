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

parser.add_argument('-Ni', '--Ninit', nargs='?', default='-1',
                    help='specify the number of initialization cycles per GCMC calculation\nDefault = min(N÷2, 10000)')

parser.add_argument('-t', '--type', nargs='?', default='grid',
                    help='specify the type of simulation you want to run on the structures\nDefault = grid')

parser.add_argument('-m', '--molecules', nargs='+', default=['xenon','krypton'],
                    help='specify the adsorbent molecules\nExample: xenon krypton CO2.\nDefault=xenon krypton')

parser.add_argument('-T', '--Temperatures', nargs='+', default=[298.0],
                    help='specify the temperatures for Raspa2 simulations \nDefault=298.0')

parser.add_argument('-p', '--pressures', nargs='+', default=['101300'],
                    help='specify the pressures in the simulation\nExample: 101300 1e6 2.5e7. Default=101300')

parser.add_argument('-Pp', '--pressureparallelism', nargs='?', default='parallel',
                    help='choose whether to let Raspa compute the different pressions sequentially in the same file ("raspa"), sequentially using each Restart file for the next batch ("sequence"), or in parallel ("parallel", default).')

parser.add_argument('-C', '--Cutoff', nargs='?', default=12,
                    help='specify the cut-off of the Raspa2 simulations\nDefault=12')

parser.add_argument('-E', '--Ewald', nargs='?', default=1e-6,
                    help='specify the Ewald precision of the Raspa2 simulations (or put 0 for no Ewald) \nDefault=1e-6')

parser.add_argument('-c', '--composition', nargs='*', default=None,
                    help='for adsorption, specify if each molecule is a cation (1) or not (0). Example: 1 0 0.\nFor coadsorption simulation, specify the composition of each adsorbent molecule. Example: 90 10. Default=None')

parser.add_argument('-r', '--radius', nargs='?', default='1.2',
                    help='specify the radius of the probe in Zeo++ calculations or the sampling sphere relative radius\nDefault: 1.2')

parser.add_argument('-rj', '--rejection', nargs='?', default='0.85',
                    help='specify the rejection condition relative radius in surface simulations \nDefault: 0.85')

parser.add_argument('-th', '--threshold', nargs='?', default='0',
                    help='reject structures with a volume above the specified threshold \nDefault=None')

parser.add_argument('-R', '--Restart', action='store_true',
                    help='specify if you want to restart from a previous state')

parser.add_argument('-sd', '--skipdone', action='store_true',
                    help='skip Raspa computations that have already been completed')

parser.add_argument('-M', '--Movie', action='store_true',
                    help='specify if you want to output the movie (for RASPA simulations)')

parser.add_argument('-o', '--output_directory', nargs='?', default='.',
                    help='specify the directory in which you want the simulation files to be installed. Default=. (current directory)')

parser.add_argument('-x', '--extra', nargs='*', type=str, default='',
                    help='extra options, given as a string to insert in Raspa2 INPUT files')

parser.add_argument('-X', '--execute', nargs='?', default='exe',
                    help='specify if you want to execute the command directly ("exe"), generate a glost list for mprun ("glost", for cluster), or prepare a slurm job ("slurm", for cluster). Default=exe')

args = parser.parse_args()

FORCE_FIELD = args.forcefield
structures_file = args.structures
option = args.type

MOLECULES = args.molecules
OUTPUT_PATH = args.output_directory
PRESSURES = args.pressures
temperatures = args.Temperatures
cutoff = float(args.Cutoff)
EwaldPrecision = float(args.Ewald)

CYCLES = int(args.Ncycles)
N_init = int(args.Ninit)

COMPOSITION = args.composition
radius = float(args.radius)
Threshold_volume = float(args.threshold)
RESTART = args.Restart
SKIPDONE = args.skipdone
PressureParallelism = args.pressureparallelism

MOVIE = args.Movie
EXTRA=' '.join(args.extra).replace(',', '\n')

nprocs = int(args.nprocesses)
ppn = int(args.procspernode)
# initialising the screening procedure
screen = Screening(structures_file, ppn, nprocs, pressures=PRESSURES, temperatures=temperatures, cutoff=cutoff,probe_radius=radius,
         force_field=FORCE_FIELD, MOLECULES=MOLECULES, type_=option, composition=COMPOSITION, cycles=CYCLES, OUTPUT_PATH=OUTPUT_PATH,
         EwaldPrecision=EwaldPrecision, Threshold_volume=Threshold_volume, RESTART=RESTART, N_init=N_init, MOVIE=MOVIE, EXTRA=EXTRA,
         SKIPDONE=SKIPDONE, PressureParallelism=PressureParallelism)

# launching the screening
if args.execute == 'glost':
    screen.glost_list()
elif args.execute == 'slurm':
    screen.slurm_job()
else:
    screen.mp_run()
