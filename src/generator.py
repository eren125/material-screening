import sys
import os
SOURCE_DIR = os.path.dirname(os.path.abspath(__file__)) 
sys.path.append(SOURCE_DIR) # So that we can import from other libraries of src/
import pandas as pd
from bash_commands import copy_dir,sed,run_shell_script
from math import ceil
# TODO add --option run

class GenerateFiles():
    def __init__(self):
        try:
            self.RASPA_DIR = os.environ['RASPA_DIR']  # read the environment variable ()
        except KeyError:
            print('RASPA_DIR is not set in your environment please run: \nsource set_environment')
        self.source_dir = os.path.dirname(os.path.abspath(__file__))
        self.current_dir = os.getcwd()
        run_shell_script("bash " + os.path.join(self.source_dir, '../data/Refresh_structures.sh'))
        self.df_mol = pd.read_csv(os.path.join(
            self.source_dir, '../data/molecules.csv'), encoding='utf-8')
        self.mol2atoms = {row['MOLECULE']: row['ATOMS']
                          for row in self.df_mol.iloc}
        self.df_structures = pd.read_csv(os.path.join(
            self.source_dir, '../data/structures.csv'), encoding='utf-8')

    def make_grid(self, FORCE_FIELD, FRAMEWORK_NAME, atom, Path_to_directory, unitcell_pattern):
        copy_dir(self.source_dir, '../Raspa_GCMC_templates/grid',
                     Path_to_directory, 'grid')
        sim_file = os.path.join(
            Path_to_directory, 'grid', 'simulation.input')
        sed("FORCE_FIELD", FORCE_FIELD, sim_file)
        sed("FRAMEWORK_NAME", FRAMEWORK_NAME, sim_file)
        sed("UnitCells 1 1 1", unitcell_pattern, sim_file)
        sed("N_ATOMS", '%d' % len(atom.split(' ')), sim_file)
        sed("ATOMS", '%s' % atom, sim_file)

        print('Generating grids for %s in %s %s' %
              (atom, FRAMEWORK_NAME, unitcell_pattern.replace('UnitCells ', '')))
        run_file = os.path.join(Path_to_directory, 'grid', 'run')
        run_shell_script(run_file)

    def generate(self, FORCE_FIELD, FRAMEWORK_NAME, MOLECULE, output_path, type_="grid", COMPOSITION=None):
        if not (FRAMEWORK_NAME+".cif" in list(self.df_structures["Structures"])):
            raise ValueError(
                "the structure you want to simulate is not in the database of the Raspa simulation directory")
        Path_to_directory = os.path.join(self.current_dir, output_path)

        # UnitCell
        path_to_unitcell = os.path.join(
            Path_to_directory, 'Unitcell', 'UnitCell.txt')
        if not os.path.exists(path_to_unitcell):
            copy_dir(self.source_dir, '../Raspa_GCMC_templates/Unitcell',
                     Path_to_directory, 'Unitcell')
            # Replace arguments accordingly
            sim_file = os.path.join(
                Path_to_directory, 'Unitcell', 'simulation.input')
            sed('FORCE_FIELD', FORCE_FIELD, sim_file)
            sed('FRAMEWORK_NAME', FRAMEWORK_NAME, sim_file)
            # Get the Unitcell size
            run_file = os.path.join(Path_to_directory, 'Unitcell', 'run')
            run_shell_script(run_file)
        # Set the UnitCell so that it is higher than the Van der Waals cutoff (12 angstrom)
        unitcell_file = open(path_to_unitcell, 'r')
        unitcell_list = unitcell_file.readlines()[0].split()
        unitcell_file.close()
        unitcell_pattern = "UnitCells %d %d %d" % tuple(
            [ceil(12/float(a)) for a in unitcell_list])

        if type_ == "ads":
            print('ADSORPTION FILES')
            # Adsorptions
            for molecule in MOLECULE:
                if not (molecule in list(self.df_mol['MOLECULE'])):
                    raise NameError('%s not in adsorbent list' % molecule)
                ATOMS = self.mol2atoms[molecule]
                ATOMS_list = ATOMS.split()
                for atom in ATOMS_list:
                    if not os.path.exists(self.RASPA_DIR + "/share/raspa/grids/%s/%s/0.100000/%s_%s_shifted.grid" % (FORCE_FIELD, FRAMEWORK_NAME, FRAMEWORK_NAME, atom)):
                        self.make_grid(FORCE_FIELD, FRAMEWORK_NAME,
                                       atom, Path_to_directory, unitcell_pattern)
                    else:
                        print('%s in %s using %s grid already calculated' %
                              (atom, FRAMEWORK_NAME, FORCE_FIELD))
                copy_dir(self.source_dir, '../Raspa_GCMC_templates/Adsorption_template',
                         Path_to_directory, 'Adsorption_'+molecule)

                sim_file = os.path.join(
                    Path_to_directory, 'Adsorption_'+molecule, 'simulation.input')
                sed("FORCE_FIELD", FORCE_FIELD, sim_file)
                sed("FRAMEWORK_NAME", FRAMEWORK_NAME, sim_file)
                sed("UnitCells 1 1 1", unitcell_pattern, sim_file)
                sed("N_ATOMS", str(len(ATOMS_list)), sim_file)
                sed("ATOMS", ATOMS, sim_file)
                sed("MOLECULE", molecule, sim_file)

        elif type_ == "grid":
            print('GRID FILES')
            ATOMS = ' '.join([self.mol2atoms[molecule]
                              for molecule in MOLECULE])
            self.make_grid(FORCE_FIELD, FRAMEWORK_NAME, ATOMS,
                           Path_to_directory, unitcell_pattern)

        elif type_ == "coad":
            print('COADSORPTION FILES')
            if not COMPOSITION:
                raise ValueError(
                    "We need composition information to run coadsorption simulation")
            if sum([int(c) for c in COMPOSITION]) != 100:
                raise ValueError("The composition doesn't amount to 100%")
            if len(MOLECULE) != 2:
                raise ValueError(
                    "We need exactly 2 molecules for a coadsorption simulation. Please change your parameters accordingly")
            # Coadsorption
            # make sure the grid / adsorption for them have been run (option to run both ?)
            composition = ' '.join(COMPOSITION)  # composition in percentage
            adsorbent = ' '.join(MOLECULE)
            grid_ATOMS = []
            for molecule in MOLECULE:
                if not (molecule in list(self.df_mol['MOLECULE'])):
                    raise NameError('%s not in adsorbent list' % molecule)
                ATOMS = self.mol2atoms[molecule]
                ATOMS_list = ATOMS.split()
                grid_ATOMS = grid_ATOMS + ATOMS_list
                for atom in ATOMS_list:
                    if not os.path.exists(self.RASPA_DIR + "/share/raspa/grids/%s/%s/0.100000/%s_%s_shifted.grid" % (FORCE_FIELD, FRAMEWORK_NAME, FRAMEWORK_NAME, atom)):
                        self.make_grid(
                            FORCE_FIELD, FRAMEWORK_NAME, atom, Path_to_directory, unitcell_pattern)
                    else:
                        print('%s in %s using %s grid already calculated' %
                              (atom, FRAMEWORK_NAME, FORCE_FIELD))
            Y = [int(c)/100 for c in COMPOSITION]

            filename = 'Coadsorption_%s_%s' % (
                adsorbent.replace(' ', '-'), composition.replace(' ', '-'))

            # Copy from template
            copy_dir(self.source_dir, '../Raspa_GCMC_templates/Coadsorption_template',
                     Path_to_directory, filename)

            # Replace arguments accordingly
            print('Replacing the arguments to fit the option given')
            sim_file = os.path.join(
                Path_to_directory, filename, 'simulation.input')
            sed('FORCE_FIELD', FORCE_FIELD, sim_file)
            sed('FRAMEWORK_NAME', FRAMEWORK_NAME, sim_file)
            sed("UnitCells 1 1 1", unitcell_pattern, sim_file)
            sed("N_ATOMS", len(grid_ATOMS), sim_file)
            sed("ATOMS", ' '.join(grid_ATOMS), sim_file)
            sed("MOLECULE_0", MOLECULE[0], sim_file)
            sed("Y_0", Y[0], sim_file)
            sed("MOLECULE_1", MOLECULE[1], sim_file)
            sed("Y_1", Y[1], sim_file)
        else:
            print("There is no such option as %s"%type_)
