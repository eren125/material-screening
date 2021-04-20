import os 
import errno
from pymatgen.io.cif import CifParser
from pymatgen.core import PeriodicSite 
import pandas as pd
import numpy as np 

class load():
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
        self.df_FF['atom_symbol'] = self.df_FF['atom type'].str.strip('_')
        self.atoms = atoms
        self.cutoff = cutoff
        self.temperature = temperature
        self.R = 8.31446261815324e-3 #kJ/K/mol
        if not os.path.exists("Energies"):
            os.mkdir("Energies")


    def evaluate_from_coordinates(self, structure_name, supercell, min_distance=2.0, energy_precision=6, energy_threshold=0.00, numpy_array=True):
        """Calculates the average energies/ minimum energy/ boltzmann average energies/ list of energies

        Args:
            structure_name (str): name of the structure as referenced in the cif file in the Raspa directory
            supercell (list): list of ints representing the size of each unicell vector for the supercell
            system of Raspa2 and Zeo++ (different from the one used in pymatgen)
        
        Returns:
        """
        if not os.path.exists("Coordinates"):
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), "Coordinates")
        coord_path = os.path.join("Coordinates",structure_name+'.csv')
        if not os.path.exists(coord_path):
            return np.nan,np.nan,np.nan

        df_voro = pd.read_csv(coord_path, sep=" ")
        df_voro = df_voro[df_voro["dist_to_nearest"]>min_distance]
        if len(df_voro)==0:
            return np.nan, np.nan, np.nan

        cif_file = self._load_cif_file(structure_name, self.RASPA_DIR)
        structure = CifParser.from_string(cif_file).get_structures(primitive=False)[0]
        structure.make_supercell(supercell)

        lattice_matrix = self._frac_to_cart_matrix(structure)
        inv_lattice_matrix = np.linalg.inv(lattice_matrix)

        df_voro["cartesian_coordinates"] = df_voro.apply(lambda row: np.array([float(row.x),float(row.y),float(row.z)]), axis=1)
        df_voro["fractional_coordinates"] = df_voro["cartesian_coordinates"].apply(lambda cart_coord: inv_lattice_matrix@np.array(cart_coord))

        accessible_mean_energy, min_energy, boltz_energy = [], [], []
        for atom_type in self.atoms:
            df_voro["pymat_site_%s"%atom_type] = df_voro["fractional_coordinates"].apply(lambda frac_coord: PeriodicSite({atom_type:1}, frac_coord, structure.lattice))
            if numpy_array:
                df_voro["Energy_%s"%atom_type] = df_voro["pymat_site_%s"%atom_type].apply(lambda site: self.lennard_jones_from_pymatgen_np(site, structure))
            else:
                df_voro["Energy_%s"%atom_type] = df_voro["pymat_site_%s"%atom_type].apply(lambda site: self.lennard_jones_from_pymatgen(site, structure))
            accessible_mean_energy.append( round(df_voro[df_voro["Energy_%s"%atom_type]<energy_threshold]["Energy_%s"%atom_type].mean(), energy_precision) )
            min_energy.append( round(df_voro["Energy_%s"%atom_type].min(), energy_precision) )
            df_voro["exp_energy_%s"%atom_type] = np.exp(-df_voro["Energy_%s"%atom_type]/(self.R*self.temperature))
            Sum_exp = df_voro["exp_energy_%s"%atom_type].sum()
            df_voro["pond_energy_%s"%atom_type] = df_voro["exp_energy_%s"%atom_type] * df_voro["Energy_%s"%atom_type] / Sum_exp
            boltz_energy.append( round(df_voro["pond_energy_%s"%atom_type].sum(), energy_precision) )
        df_voro.to_csv('Energies/%s.csv'%structure_name)
        return accessible_mean_energy, min_energy, boltz_energy


    def evaluate_from_surface(self, structure_name, N_sample = 500, min_distance=2.0, energy_precision=6, energy_threshold=0.00, numpy_array=True):
        """A function to sample the accessible surfaces and calculate the energies associated
        Args:
            structure_name (str): the name of the structure (need to be present in the RASPA directory)
            
        """
        cif_file = self._load_cif_file(structure_name, self.RASPA_DIR)
        structure = CifParser.from_string(cif_file).get_structures(primitive=False)[0]
        N_atoms = len(structure)

        lattice_matrix = structure.lattice.matrix

        accessible_mean_energy, min_energy, boltz_energy = [], [], []

        for atom_type_g in self.atoms:
            E_list = []

            row_g = self.df_FF[self.df_FF['atom type'] == atom_type_g]
            sigma_g = row_g['sigma'].iloc[0]
            epsilon_g = row_g['epsilon'].iloc[0]

            center_indices, neighbor_indices, images, distances = structure.get_neighbor_list(self.cutoff + sigma_g, exclude_self=False)

            for index in range(N_atoms):
                atom_h = structure[index]
                atom_type_h = atom_h.specie.symbol
                row_h = self.df_FF[self.df_FF['atom type'].str.strip('_') == atom_type_h]
                sigma = (sigma_g + row_h['sigma'].iloc[0]) / 2
                radius = sigma * pow(2,1/6)

                index_slice = np.where(center_indices==index)[0]
                deb, fin = index_slice[0], index_slice[-1]

                neighbors_i = neighbor_indices[deb:fin]
                images_i = images[deb:fin]

                cart_coord_neighbors = np.array([np.dot([structure[j].a + im[0], structure[j].b + im[1], structure[j].c + im[2]], lattice_matrix) for j, im in zip(neighbors_i, images_i)])

                atom_type_neighbors = np.array([structure[j].specie.symbol for j in neighbors_i])

                df_temp = pd.DataFrame(data={'atom_symbol': atom_type_neighbors})
                df_temp = pd.merge(df_temp, self.df_FF, how='left', on='atom_symbol')

                epsilon = np.sqrt( epsilon_g * df_temp['epsilon'].to_numpy() )
                sigma = ( sigma_g + df_temp['sigma'].to_numpy() ) /2
                # TODO use an array instead of this loop
                for sample in range(N_sample):
                    u, v, w = np.random.normal(0,1), np.random.normal(0,1), np.random.normal(0,1)
                    d = np.sqrt(u**2 + v**2 + w**2)
                    cart_coord_adsorbent = atom_h.coords + radius*np.array([u,v,w])/d
                    dist_absorbent_neighbor = self._distance(cart_coord_adsorbent, cart_coord_neighbors)
                    # cutoff
                    dist_cutoff = np.extract(dist_absorbent_neighbor<self.cutoff, dist_absorbent_neighbor)
                    epsilon_cutoff = np.extract(dist_absorbent_neighbor<self.cutoff, epsilon)
                    sigma_cutoff = np.extract(dist_absorbent_neighbor<self.cutoff, sigma)

                    E = self._calculate_lj(epsilon_cutoff, sigma_cutoff, dist_cutoff)
                    shift = self._calculate_lj(epsilon_cutoff, sigma_cutoff, self.cutoff)

                    LJ_energy = self.R * (E-shift)
                    E_list.append(LJ_energy)
            E_list_np = np.array(E_list)
            E_exp = np.exp(-E_list_np/(self.R*self.temperature))
            S = np.sum(E_exp)
            E_Boltz = E_exp * E_list_np / S
            accessible_energy = np.extract(E_list_np<energy_threshold,E_list_np)
            if len(accessible_energy)>0:
                accessible_mean_energy.append(round(np.mean(accessible_energy), energy_precision))
                min_energy.append(round(np.min(accessible_energy), energy_precision))
            else:
                accessible_mean_energy.append(np.nan)
                min_energy.append(round(np.min(E_list_np), energy_precision))
            boltz_energy.append(round(np.sum(E_Boltz), energy_precision))
        return accessible_mean_energy, min_energy, boltz_energy


    def lennard_jones_from_pymatgen_np(self, atom_g, structure_h, shifted=True):
        """Calculates the VdW interaction using a LJ potential

        Args:
            atom_g (pymatgen.core.PeriodicSite): element type of the guest atom
            structure_h (pymatgen.core.Structure): element type of the host atom
            distance (float): distance between the atoms
            df_FF (pd.DF): a pandas DataFrame giving the LJ parameters of the corresponding force field

        Returns:
            Lennard Jones Interaction Energy in K according to the LJ parameters of df_FF
        """

        shift = 0
        row_g = self.df_FF[self.df_FF['atom type'] == atom_g.specie.symbol]
        epsilon_g = row_g['epsilon'].iloc[0]
        sigma_g = row_g['sigma'].iloc[0]
        # Improve using from pymatgen.optimization.neighbors import find_points_in_spheres
        neighbors = structure_h.get_neighbors(atom_g, self.cutoff)

        distance = np.array([atom_g.distance(atom_h) for atom_h in neighbors])
        df_temp = pd.DataFrame(data={"atom_symbol":[atom_h.specie.symbol for atom_h in neighbors]})
        df_temp = pd.merge(df_temp, self.df_FF, how='left', on='atom_symbol')

        epsilon = np.sqrt(epsilon_g * df_temp['epsilon'].to_numpy())
        sigma = (sigma_g + df_temp['sigma'].to_numpy()) / 2

        E = self._calculate_lj(epsilon, sigma, distance)
        if shifted:
            shift = self._calculate_lj(epsilon, sigma, self.cutoff)
        return self.R * (E - shift)


    def lennard_jones_from_pymatgen(self, atom_g, structure_h, shifted=True):
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
        neighbors = structure_h.get_neighbors(atom_g, self.cutoff)
        for atom_h in neighbors:
            distance = atom_h.distance(atom_g)
            row_h = self.df_FF[self.df_FF['atom type'].str.strip('_') == atom_h.specie.symbol]
            epsilon = np.sqrt(row_g['epsilon'].iloc[0] * row_h['epsilon'].iloc[0])
            sigma = (row_g['sigma'].iloc[0] + row_h['sigma'].iloc[0]) / 2
            E += self._calculate_lj(epsilon, sigma, distance)
            if shifted:
                shift += self._calculate_lj(epsilon, sigma, self.cutoff)
        return self.R * (E - shift)


    @staticmethod
    def _calculate_lj(epsilon, sigma, distance):
        """Function to calculate the Lennard-Jones interactions using epsilon, sigma and the distance
        """
        return np.sum( 4 * epsilon * ( (sigma / distance)**12 - (sigma / distance)**6 ) )

    @staticmethod
    def _load_cif_file(structure_name, RASPA_DIR):
        """ Load the cif file corresponding to structure_name in the RASPA directory specified by RASPA_DIR
        """
        cif_file = open(os.path.join(RASPA_DIR,'share/raspa/structures/cif',structure_name+".cif"),"r").read()

        if '_atom_site_label' in cif_file:
            if '_atom_site_type_symbol' in cif_file:
                cif_file = cif_file.replace('_atom_site_label','_atom_site_lab')
                cif_file = cif_file.replace('_atom_site_type_symbol','_atom_site_label')
        else:
            if '_atom_site_type_symbol' in cif_file:
                cif_file = cif_file.replace('_atom_site_type_symbol','_atom_site_label')
            else:
                raise KeyError("_atom_site_type_symbol and _atom_site_label not in the cif file, please check the format")
        return cif_file

    @staticmethod
    def _frac_to_cart_matrix(structure):
        """ Using the crystallography text book formula to convert fractionnal
        coordinates to cartesian coordinates
        Args:
            structure (pymatgen.core.Structure): pymatgen parsed structure extracted from cif file. It contains the structures' lattice properties used to determine the final lattice_matrix
        Output: 
            lattice_matrix
        """
        a, b, c = structure.lattice.abc
        alpha, beta, gamma = structure.lattice.angles
        alpha = np.pi*alpha/180
        beta = np.pi*beta/180
        gamma = np.pi*gamma/180
        n = ( np.cos(alpha) - np.cos(gamma)*np.cos(beta) ) / np.sin(gamma)
        return np.array([[a, round(b*np.cos(gamma),12), round(c*np.cos(beta),12)],
                         [0.0, round(b*np.sin(gamma),12), round(c*n,12)],
                         [0.0, 0.0, round(c*np.sqrt(np.sin(beta)**2 - n**2),12)]
                        ])

    @staticmethod
    def _distance(coord_1, coord_array_2):
        """ Calculates the distances between a point of coordinate coord_1 and all the points of coord_array_2
        """
        return np.sqrt( (coord_1[0]-coord_array_2[:,0])**2 + (coord_1[1]-coord_array_2[:,1])**2 + (coord_1[2]-coord_array_2[:,2])**2 )
