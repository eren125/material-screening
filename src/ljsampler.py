import os 
from pymatgen.io.cif import CifParser
from pymatgen.core import PeriodicSite 
import pandas as pd
import numpy as np 

class load_coordinates():
    def __init__(self, forcefield = "UFF"):
        """Initialise the forcefield and creates a dataframe of forcefield used throughout the code
        """

        self.RASPA_DIR = os.environ['RASPA_DIR']
        FF_filelines = open(os.path.join(self.RASPA_DIR,'share/raspa/forcefield',forcefield,"force_field_mixing_rules.def"),"r").readlines()
        columns = FF_filelines[6].strip('#').strip().split(', ')
        FF_content = [l.split() for l in FF_filelines[7:-2]]
        self.df_FF = pd.DataFrame(FF_content, columns=columns).astype({'epsilon':float, "sigma":float})
        if self.df_FF.loc[0,'interaction form'] != 'lennard-jones':
            raise ValueError("Unsupported interation form: only Lennard-Jones interactions are supported")
        

    def sample(self, structure_name, unitcell, box_matrix):
        """Calculates the average energies/ minimum energy/ boltzmann average energies/ list of energies (print in csv??)
        """

        df_voro = pd.read_csv(os.path.join("Coordinates",structure_name+'.csv'), sep=" ")
        cif_file = open(os.path.join(self.RASPA_DIR,'share/raspa/structures/cif',structure_name+".cif"),"r").read()
        cif_file = cif_file.replace('_atom_site_type_symbol','_atom_site_label')

        pass

    def lennard_jones(self, atom_g, atom_h, distance, cutoff, shifted=True):
        """Calculates the VdW interaction using a LJ potential
        Args:
            atom_g (str): element type of the guest atom
            atom_h (str): element type of the host atom
            distance (float): distance between the atoms
            df_FF (pd.DF): a pandas DataFrame giving the LJ parameters of the corresponding force field
        Returns:
            Lennard Jones Interaction Energy in K according to the LJ parameters of df_FF
        """

        if 0 < distance < cutoff:
            row_g = self.df_FF[self.df_FF['atom type'] == atom_g]
            row_h = self.df_FF[self.df_FF['atom type'].str.strip('_') == atom_h]
            epsilon = np.sqrt(row_g['epsilon'].iloc[0] * row_h['epsilon'].iloc[0])
            sigma = (row_g['sigma'].iloc[0] + row_h['sigma'].iloc[0]) / 2
            E = 4 * epsilon * ( pow(sigma,12)/pow(distance,12) - pow(sigma,6)/pow(distance,6) )
            if shifted:
                return E - 4 * epsilon * ( pow(sigma,12)/pow(cutoff,12) - pow(sigma,6)/pow(cutoff,6) )
            else:
                return E
        elif distance <= 0:
            raise ValueError("Negative distance: please chack the data fed in this function")
        else:
            return 0.0
