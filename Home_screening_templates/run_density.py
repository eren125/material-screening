import sys
import os
import ase
import ase.io


def calculate_density(cell):

    vol = cell.get_volume()
    mass = sum(cell.get_masses())

    # Return in g.cm^-3
    return mass / vol * 1e30 / 6.02214076e+23 / 1e6


##############################################
# Main program follows

if len(sys.argv) < 1:
    prog = os.path.basename(sys.argv[0])
    print(f'Usage: {prog} file.cif')
    sys.exit(1)

try:
    cell = ase.io.read(sys.argv[1])
except Exception as e:
    print('Could not read CIF structure from file: ' + sys.argv[1])
    print('Error message is: ' + str(e))
    sys.exit(1)

rho = calculate_density(cell)
os.system('echo %s,%s >> density'%(sys.argv[1],rho))

sys.exit(0)
