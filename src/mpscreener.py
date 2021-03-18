import sys
import os
SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SOURCE_DIR)

import pandas as pd
import numpy as np
import multiprocessing as mp
from time import time  # Temporary for testing
from bash_commands import copy_dir, copy_file, sed, bash


# TODO
# catch errors > type errors ?
# ADD cutoff as a variable for raspa simulations


class Screening():
    def __init__(self, structures_file, FORCE_FIELD, MOLECULES, PRESSURES, nprocs, OUTPUT_PATH=".", TEMP=298.0, 
    probe_radius=1.2, Threshold_volume=20, procs_per_node=48, type_='grid', COMPOSITION=None, CYCLES=2000,setup=True):
        """A class for screening purposes using Raspa2 for molecular simulations
        (the path of the Raspa directory have to be exported as en environement variable)

        Args:
            structures_file    (str): relative path to the csv file containing all the structures, 
                                      the headers must contain "Structures"
            FORCE_FIELD        (str): force field used for the molecular simulations, 
                                      it must be defined in the Raspa directory 
            MOLECULES         (list): list of molecules (str) to be adsorbed on the materials
            PRESSURES         (list): list of pressures (float) in pascal to be simulated
            nprocs             (int): number of processes to be run at the same time
            OUTPUT_PATH        (str): relative path to where the outputs of the simulations will be printed
            TEMP             (float): temperature in kelvin for the Raspa2 simulations (default = 298.0 K)
            probe_radius     (float): radius of the probe in angström considered in Zeo++ simulations
            Threshold_volume (float): threshold at which the volume is considered too big for grid calculations
                                      it consumes too much RAM Memory for the machine currently used
            procs_per_node     (int): number of processors available in each node for a given calculator
            type_              (str): type of simulation the user wants to carry out. Currently available:
                                      grid calculation via 'grid', GCMC via 'ads' 'coad', NVT MC via 'ent',
                                      Widom's insertion via 'widom', helium void fraction via 'vf', Zeo++ via 'surface' 'volume' 'pore' 'channel' and global information via 'info'.
            COMPOSITION       (list): list of mole fractions for the molecules defined in MOLECULES
                                      must have the same length as MOLECULES, and sum must equal to 1
            CYCLES             (int): number of prodution cycles used for the Raspa2 simulations
                                      equilibration cycles are fixed to 10k for now (to improve)
            setup             (bool): boolean variable to determine whether to copy the input files 
                                      (obsolete in next version: the input will be parsed)

        self variables:
            NODES              (str): names of the nodes allocated by the supercalculator
                                      if you are running on a local computer 'NODES' does not exist
                                      it will then be set to an empty string by default
            df_mol           (pd.df): pandas dataframe containing the available molecules and 
                                      the atoms it is composed of
            mol2atoms         (dict): dictionnary with molecule names as keys and list of atoms associated as values
            df_structures    (pd.df): dataframe containing the structures to be simulated in this run

        self methods:
            run   : function that takes the framework's name and the smallest unitcell for a 12 angstrom cut-off
                    and runs the corresponding simulation according to the run bash file 
                    (TODO to be improved using wrapsa2)
            run_mp: t

        """
        try:
            self.NODES = os.environ['NODES']
        except:
            self.NODES = ''
        self.NODES = self.NODES.split()
        if procs_per_node > nprocs:
            raise ValueError('More processes at a time than proccessors available!!!')
        self.nprocs = nprocs
        available_types = ['ads', 'ent', 'widom', 'coad', 'grid', 'info','vf','surface','volume','pore','channel']
        if type_ not in available_types:
            raise ValueError(('%s not an option yet. Please choose between: ' +
                ', '.join(['%s']*len(available_types))) % tuple([type_]+available_types))
        self.type_ = type_

        # conversion table molecule to atoms and vice versa
        self.df_mol = pd.read_csv(os.path.join(SOURCE_DIR, '../data/molecules.csv'), encoding='utf-8')
        self.mol2atoms = {row['MOLECULE']: row['ATOMS'].split() for index,row in self.df_mol.iterrows()}
        if not all([M in self.df_mol['MOLECULE'] for M in MOLECULES]):
            raise ValueError("The molecules mentioned are not supported by the code, see the data directory in the source directory")
        self.MOLECULES = MOLECULES

        # Structure dataframe construction
        current_dir = os.getcwd()
        self.df_structures = pd.read_csv(os.path.join(current_dir, structures_file), encoding='utf-8')
        self.df_structures = self.df_structures[['Structures']]
        self.df_structures['STRUCTURE_NAME'] = self.df_structures['Structures'].str.replace('.cif','')
        if type_ == 'info':
            self.data = self.df_structures[['STRUCTURE_NAME','Structures']].to_records(index=False)
        else:
            # import information on all available structures
            df_info = pd.read_csv(os.path.join(SOURCE_DIR, "../data/info.csv"), encoding='utf-8')
            df = pd.merge(self.df_structures[['STRUCTURE_NAME']], df_info[['STRUCTURE_NAME','UnitCell','Volume [nm^3]']],how="inner", on="STRUCTURE_NAME")
            # remove small volumes
            df = df[df['Volume [nm^3]'] <= Threshold_volume]

        # Remove this part or use it to chose between the raspa functions defined in raspa2.py
        Path_to_directory = os.path.join(current_dir, OUTPUT_PATH)
        if setup==True:
            if type_ in ['grid','ads','coad','ent','widom','vf']: 
                path_to_Scripts = os.path.join(Path_to_directory, 'Scripts')
                if not os.path.exists(path_to_Scripts):
                    os.mkdir(path_to_Scripts)
                copy_file(SOURCE_DIR,"../Raspa_screening_templates/INPUT_%s"%type_,Path_to_directory,"INPUT")
                copy_file(SOURCE_DIR,"../Raspa_screening_templates/run",Path_to_directory,"run")
                if type_ != 'grid':
                    copy_file(SOURCE_DIR,"../Raspa_screening_templates/data_%s.sh"%type_,Path_to_directory,"data.sh")
            elif type_ == "info":
                path_to_Scripts = os.path.join(Path_to_directory, 'Scripts')
                if not os.path.exists(path_to_Scripts):
                    os.mkdir(path_to_Scripts)
                copy_file(SOURCE_DIR,"../Raspa_screening_templates/INPUT_%s"%type_,Path_to_directory,"INPUT")
                copy_file(SOURCE_DIR,"../Raspa_screening_templates/run_%s"%type_,Path_to_directory,"run")
            else :
                path_to_output = os.path.join(Path_to_directory, 'Output')
                if not os.path.exists(path_to_output):
                    os.mkdir(path_to_output)
                copy_file(SOURCE_DIR,"../Zeo++_screening_templates/run_%s"%type_,Path_to_directory,"run")
                df['UnitCell'] = probe_radius
                
            self.data = df[['STRUCTURE_NAME','UnitCell']].to_records(index=False)

            self.path_to_run = os.path.join(Path_to_directory,"run")
            path_to_data = os.path.join(Path_to_directory,"data.sh")
            # sed
            path_to_INPUT = os.path.join(Path_to_directory, 'INPUT')
            sed("FORCE_FIELD",FORCE_FIELD,path_to_INPUT)
            if CYCLES:
                sed("N_cycles",CYCLES,path_to_INPUT)
                sed("NCYCLES",CYCLES,self.path_to_run)

            ATOMS = []
            for molecule in MOLECULES:
                if not (molecule in list(self.df_mol['MOLECULE'])):
                    raise NameError('%s not in adsorbent list' % molecule)
                atoms = self.mol2atoms[molecule]
                ATOMS += atoms

            sed("N_ATOMS", len(ATOMS), path_to_INPUT)
            sed("ATOMS",' '.join(ATOMS),path_to_INPUT)
            sed("PRESSURES",' '.join(PRESSURES),path_to_INPUT)
            if len(MOLECULES) == 1:
                sed("MOLECULE",MOLECULES[0],path_to_INPUT)
                sed("MOLECULE",MOLECULES[0],path_to_data)
            elif len(MOLECULES) == 2:
                if COMPOSITION:
                    Y = [int(c)/100 for c in COMPOSITION]
                    if sum(Y) != 1.0:
                        raise ValueError("Composition does not add up to 100%%")
                    sed("MOLECULE_1",MOLECULES[0],path_to_INPUT)
                    sed("Y_1",Y[0],path_to_INPUT)
                    sed("MOLECULE_2",MOLECULES[1],path_to_INPUT)
                    sed("Y_2",Y[1],path_to_INPUT)

    def run(self, inputs):
        if len(self.NODES) == 0:
            command = "bash %s %s \\\"%s\\\""%(self.path_to_run,FRAMEWORK_NAME,UNITCELL)
        else:
            worker = int(mp.current_process()._identity[0])
            nnode = len(self.NODES)
            index = (worker-1)%nnode
            FRAMEWORK_NAME, UNITCELL = inputs
            command = "ssh %s \"bash %s %s \\\"%s\\\"\""%(self.NODES[index],self.path_to_run,FRAMEWORK_NAME,UNITCELL)
        return os.system(command)

    def mp_run(self):
        t0 = time()
        print("%s simulation of %s"%(self.type_,' '.join(self.MOLECULES)))
        with mp.Pool(processes=self.nprocs) as p:
            p.map(self.run, [(FRAMEWORK_NAME,UNITCELL) for FRAMEWORK_NAME,UNITCELL in self.data])
        print("SIMULATIONS COMPLETED after %.2f hour(s)"%((time()-t0)/3600))
