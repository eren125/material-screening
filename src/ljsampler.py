import os 
from pymatgen.io.cif import CifParser
from pymatgen.core import PeriodicSite 
import pandas as pd
import numpy as np 

class load_coordinates():
    def __init__(self, atoms, forcefield = "UFF", temperature=298.0, cutoff=12.0):
        """Initialise the forcefield and creates a dataframe of forcefield used throughout the code
        Args:
            atoms (list): a list of strings precising the adsorbent atoms to use in the energy calculation
            forcefield (str): the force field name used for the energy calculation (need to be specified in the raspa directory)
            temperature (float): the temperature used for the simulation
            cutoff (float): the cutoff used in the energy calculation (12 angstrom is the default)

        """

        self.RASPA_DIR = os.environ['RASPA_DIR']
        FF_filelines = open(os.path.join(self.RASPA_DIR,'share/raspa/forcefield',forcefield,"force_field_mixing_rules.def"),"r").readlines()
        columns = FF_filelines[6].strip('#').strip().split(', ')
        FF_content = [l.split() for l in FF_filelines[7:-2]]
        self.df_FF = pd.DataFrame(FF_content, columns=columns).astype({'epsilon':float, "sigma":float})
        if self.df_FF.loc[0,'interaction form'] != 'lennard-jones':
            raise ValueError("Unsupported interation form: only Lennard-Jones interactions are supported")
        self.atoms = atoms
        self.cutoff = cutoff
        self.temperature = temperature
        self.R = 8.31446261815324e-3 #kJ/K/mol
        if not os.path.exists("Energies"):
            os.mkdir("Energies")

    def sample(self, structure_name, unitcell, lattice_matrix, min_distance=2.0):
        """Calculates the average energies/ minimum energy/ boltzmann average energies/ list of energies

        Args:
            structure_name (str): name of the structure as referenced in the cif file in the Raspa directory
            unitcell (str): supercell unitcell to prevent self interaction at 12 angstrom cut off
            lattice_matrix (numpy.2Darray): 3x3 array representing the unitcell vectors in the cartesian coordinate 
            system of Raspa2 and Zeo++ (different from the one used in pymatgen)

        """

        df_voro = pd.read_csv(os.path.join("Coordinates",structure_name+'.csv'), sep=" ")
        df_voro = df_voro[df_voro["dist_to_nearest"]>min_distance]
        cif_file = open(os.path.join(self.RASPA_DIR,'share/raspa/structures/cif',structure_name+".cif"),"r").read()
        cif_file = cif_file.replace('_atom_site_type_symbol','_atom_site_label')
        parser = CifParser.from_string(cif_file)
        structure = parser.get_structures()[0]
        structure.make_supercell([int(u) for u in unitcell.split()])
        df_voro["cartesian_coordinates"] = df_voro.apply(lambda row: np.array([float(row.x),float(row.y),float(row.z)]), axis=1)
        inv_lattice_matrix = np.linalg.inv(lattice_matrix)
        df_voro["fractional_coordinates"] = df_voro["cartesian_coordinates"].apply(lambda coord: np.matmul(inv_lattice_matrix, np.array(coord)))

        accessible_mean_energy, min_energy, boltz_energy = [], [], []
        for molecule in self.atoms:
            df_voro["pymat_site_%s"%molecule] = df_voro["fractional_coordinates"].apply(lambda frac_coord: PeriodicSite({molecule:1}, frac_coord, structure.lattice))
            df_voro["Energy_%s"%molecule] = df_voro["pymat_site_%s"%molecule].apply(lambda site: self.lennard_jones(site, structure, self.cutoff))
            accessible_mean_energy.append(df_voro[df_voro["Energy_%s"%molecule]<0]["Energy_%s"%molecule].mean())
            min_energy.append(df_voro["Energy_%s"%molecule].min())
            df_voro["exp_energy_%s"%molecule] = np.exp(-df_voro["Energy_%s"%molecule]/(self.R*self.temperature))
            Sum_exp = df_voro["exp_energy_%s"%molecule].sum()
            df_voro["pond_energy_%s"%molecule] = df_voro["exp_energy_%s"%molecule] * df_voro["Energy_%s"%molecule] / Sum_exp
            boltz_energy.append(df_voro["pond_energy_%s"%molecule].sum())
        df_voro.to_csv('Energies/%s.csv'%structure_name)
        return accessible_mean_energy, min_energy, boltz_energy


    def lennard_jones(self, atom_g, structure_h, cutoff, shifted=True):
        """Calculates the VdW interaction using a LJ potential

        Args:
            atom_g (pymatgen.core.PeriodicSite): element type of the guest atom
            structure_h (pymatgen.core.Structure): element type of the host atom
            distance (float): distance between the atoms
            df_FF (pd.DF): a pandas DataFrame giving the LJ parameters of the corresponding force field

        Returns:
            Lennard Jones Interaction Energy in K according to the LJ parameters of df_FF

        """

        E = 0.0
        shift = 0.0
        row_g = self.df_FF[self.df_FF['atom type'] == atom_g.specie.symbol]
        for atom_h in structure_h:
            distance = atom_h.distance(atom_g)
            if  distance < cutoff:
                row_h = self.df_FF[self.df_FF['atom type'].str.strip('_') == atom_h.specie.symbol]
                epsilon = np.sqrt(row_g['epsilon'].iloc[0] * row_h['epsilon'].iloc[0])
                sigma = (row_g['sigma'].iloc[0] + row_h['sigma'].iloc[0]) / 2
                E += 4 * epsilon * ( (pow(sigma,12)/pow(distance,12)) - (pow(sigma,6)/pow(distance,6)) )
                if shifted:
                    shift += 4 * epsilon * ( (pow(sigma,12)/pow(cutoff,12)) - (pow(sigma,6)/pow(cutoff,6)) )
        if shifted:
            return self.R * (E - shift)
        else:
            return self.R * E
