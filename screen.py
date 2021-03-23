#!/home/emmanuel/.venv/bin/python
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

parser.add_argument('-n', '--nprocesses', nargs='?', default='24',
                    help='specify the number of processes to run at a time\nDefault = 24')

parser.add_argument('-ppn', '--procspernode', nargs='?', default='48',
                    help='specify the number of processors per node in the cluster\nDefault = 48')

parser.add_argument('-N', '--Ncycles', nargs='?', default='1000',
                    help='specify the number of cycles per GCMC calculation\nDefault = 1000')

parser.add_argument('-t', '--type', nargs='?', default='grid',
                    help='specify the type of simulation you want to run on the structures\nDefault = grid')

parser.add_argument('-r', '--radius', nargs='?', default='1.2',
                    help='specify the radius of the probe in Zeo++ calculations\nDefault: 1.2 angstrom')

parser.add_argument('-m', '--molecules', nargs='+', default=['xenon','krypton'],
                    help='specify the adsorbent molecules\nExample: xenon krypton CO2.\nDefault=xenon krypton')

parser.add_argument('-p', '--pressures', nargs='+', default=['101300'],
                    help='specify the pressures in the simulation\nExample: xenon krypton CO2. Default=101300')

parser.add_argument('-c', '--composition', nargs=2, default=None,
                    help='specify the composition of each adsorbent molecule for coadsorption simulation\nExample: 90 10. Default=None')

parser.add_argument('-pos', '--positions', nargs='?', default=None,
                    help='specify a csv file with the coordinates of the adsorbents for single point simulations.\n Default=None')

parser.add_argument('-o', '--output_directory', nargs='?', default='.',
                    help='specify the directory in which you want the simulation files to be installed. Default=. (current directory)')

args = parser.parse_args()

FORCE_FIELD = args.forcefield
structures_file = args.structures
option = args.type

MOLECULES = args.molecules
OUTPUT_PATH = args.output_directory
PRESSURES = args.pressures

CYCLES = int(args.Ncycles)
COMPOSITION = args.composition
positions = args.positions
radius = float(args.radius)

nprocs = int(args.nprocesses)
ppn = int(args.procspernode)
screen = Screening(structures_file, FORCE_FIELD, MOLECULES, nprocs, pressures=PRESSURES, OUTPUT_PATH=OUTPUT_PATH,probe_radius=radius, 
         procs_per_node=ppn, type_=option, composition=COMPOSITION, cycles=CYCLES, positions=positions)

screen.mp_run()
