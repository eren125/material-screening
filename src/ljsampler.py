import os 
import errno
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
        if not os.path.exists("Coordinates"):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), "Coordinates")
        self.atoms = atoms
        self.cutoff = cutoff
        self.temperature = temperature
        self.R = 8.31446261815324e-3 #kJ/K/mol
        if not os.path.exists("Energies"):
            os.mkdir("Energies")

    def evaluate(self, structure_name, supercell, min_distance=2.7, energy_precision=2, energy_threshold=0.00):
        """Calculates the average energies/ minimum energy/ boltzmann average energies/ list of energies

        Args:
            structure_name (str): name of the structure as referenced in the cif file in the Raspa directory
            supercell (list): list of ints representing the size of each unicell vector for the supercell
            system of Raspa2 and Zeo++ (different from the one used in pymatgen)
        
        Returns:


        """

        coord_path = os.path.join("Coordinates",structure_name+'.csv')
        if os.path.exists(coord_path):
            df_voro = pd.read_csv(coord_path, sep=" ")
            df_voro = df_voro[df_voro["dist_to_nearest"]>min_distance]
            if len(df_voro)==0:
                return np.nan, np.nan, np.nan
            cif_file = open(os.path.join(self.RASPA_DIR,'share/raspa/structures/cif',structure_name+".cif"),"r").read()

            if '_atom_site_label' in cif_file:
                if '_atom_site_type_symbol' in cif_file:
                    cif_file = cif_file.replace('_atom_site_label','_atom_site_lab')
                    cif_file = cif_file.replace('_atom_site_type_symbol','_atom_site_label')
            else:
                if '_atom_site_type_symbol' in cif_file:
                    cif_file = cif_file.replace('_atom_site_type_symbol','_atom_site_label')
                else:
                    raise KeyError("_atom_site_type_symbol and _atom_site_label not in the cif file, please check the format")

            parser = CifParser.from_string(cif_file)
            structure = parser.get_structures(primitive=False)[0]
            structure.make_supercell(supercell)

            a, b, c = structure.lattice.abc
            alpha, beta, gamma = structure.lattice.angles
            alpha = np.pi*alpha/180
            beta = np.pi*beta/180
            gamma = np.pi*gamma/180
            n = ( np.cos(alpha) - np.cos(gamma)*np.cos(beta) ) / np.sin(gamma)
            lattice_matrix = np.array([
                [a, round(b*np.cos(gamma),12), round(c*np.cos(beta),12)],
                [0.0, round(b*np.sin(gamma),12), round(c*n,12)],
                [0.0, 0.0, round(c*np.sqrt(np.sin(beta)**2 - n**2),12)] ])
            inv_lattice_matrix = np.linalg.inv(lattice_matrix)

            df_voro["cartesian_coordinates"] = df_voro.apply(lambda row: np.array([float(row.x),float(row.y),float(row.z)]), axis=1)
            df_voro["fractional_coordinates"] = df_voro["cartesian_coordinates"].apply(lambda coord: np.matmul(inv_lattice_matrix, np.array(coord)))

            accessible_mean_energy, min_energy, boltz_energy = [], [], []
            for molecule in self.atoms:
                df_voro["pymat_site_%s"%molecule] = df_voro["fractional_coordinates"].apply(lambda frac_coord: PeriodicSite({molecule:1}, frac_coord, structure.lattice))
                df_voro["Energy_%s"%molecule] = df_voro["pymat_site_%s"%molecule].apply(lambda site: self.lennard_jones(site, structure, self.cutoff))
                accessible_mean_energy.append( round(df_voro[df_voro["Energy_%s"%molecule]<energy_threshold]["Energy_%s"%molecule].mean(), energy_precision) )
                min_energy.append( round(df_voro["Energy_%s"%molecule].min(), energy_precision) )
                df_voro["exp_energy_%s"%molecule] = np.exp(-df_voro["Energy_%s"%molecule]/(self.R*self.temperature))
                Sum_exp = df_voro["exp_energy_%s"%molecule].sum()
                df_voro["pond_energy_%s"%molecule] = df_voro["exp_energy_%s"%molecule] * df_voro["Energy_%s"%molecule] / Sum_exp
                boltz_energy.append( round(df_voro["pond_energy_%s"%molecule].sum(), energy_precision) )
            df_voro.to_csv('Energies/%s.csv'%structure_name)
            return accessible_mean_energy, min_energy, boltz_energy
        else:
            return np.nan,np.nan,np.nan


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
        neighbors = structure_h.get_neighbors(atom_g, cutoff)
        for atom_h in neighbors:
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
