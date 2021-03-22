import sys
import os
SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SOURCE_DIR)

from time import time
import numpy as np
import pandas as pd

import multiprocessing as mp

from bash_commands import copy_dir, copy_file, sed, bash

SIMULATION_TYPES = {"RASPA2" : ['grid', 'ads', 'coad', 'ent', 'widom', 'vf', 'point'],
                    "INFO"   : ['info'],
                    "ZEO++"  : ['surface', 'volume', 'pore', 'channel', 'voronoi']
                   }


# TODO
# Add single point simulation logic
# Chain several sp sim
# Move procs and procs_per_node to mp_run?
# ADD cutoff as a variable for raspa simulations

class Screening():
    def __init__(self, structures_file, force_field, MOLECULES, PRESSURES, nprocs, OUTPUT_PATH=".", TEMP=298.0, 
    probe_radius=1.2, Threshold_volume=20, procs_per_node=48, type_='grid', COMPOSITION=None, positions=None, cycles=2000, setup=True):
        """A class for screening purposes using Raspa2 for molecular simulations
        (the path of the Raspa directory have to be exported as en environement variable)

        Args:
            structures_file    (str): relative path to the csv file containing all the structures, 
                                      the headers must contain "Structures"
            force_field        (str): force field used for the molecular simulations, 
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
            cycles             (int): number of prodution cycles used for the Raspa2 simulations
                                      equilibration cycles are fixed to 10k for now (to improve)
            setup             (bool): boolean variable to determine whether to copy the input files 
                                      (obsolete in next version: the input will be parsed)

        self variables:
            NODES              (str): names of the nodes allocated by the supercalculator
                                      if you are running on a local computer 'NODES' does not exist
                                      it will then be set to an empty string by default
            mol2atoms         (dict): dictionnary with molecule names as keys and list of atoms associated as values
            data

        self methods:
            run   : function that takes the framework's name and the smallest unitcell for a 12 angstrom cut-off
                    and runs the corresponding simulation according to the run bash file 
            run_mp: function that calls `run` multiple times in parallel. mp.Pool distributes the jobs to $nprocs workers
                    so that every job  
        """

        ### Initialisation of class objects and error catching ###
        try:
            self.NODES = os.environ['NODES']
        except:
            self.NODES = ''
        self.NODES = self.NODES.split()
        if len(self.NODES)==0:
            if procs_per_node > nprocs:
                raise ValueError('More processes at a time than proccessors available!!!')
        else:
            if len(self.NODES)*procs_per_node > nprocs:
                raise ValueError('More processes at a time than proccessors available!!!')
        self.nprocs = nprocs
        available_types = ['ads', 'ent', 'widom', 'coad', 'grid', 'info','vf','surface','volume','pore','channel']

        if type_ not in available_types:
            raise ValueError(('%s not an option yet. Please choose between: ' +
                ', '.join(['%s']*len(available_types))) % tuple([type_]+available_types))

        if (type_=='point') and (not positions):
            raise ValueError("Please specify the path to a csv file with positions. See examples in data/positions_sample.csv.")

        if (type_=='coad') and (not COMPOSITION):
            raise ValueError("Please specify a composition, for example: 20 80")

        if COMPOSITION:
            if sum([float(c) for c in COMPOSITION]) != 100.0:
                raise ValueError("Composition does not add up to 100%%")
        self.type_ = type_

        df_mol = pd.read_csv(os.path.join(SOURCE_DIR, '../data/molecules.csv'), encoding='utf-8')
        if not all([molecule in list(df_mol['MOLECULE']) for molecule in MOLECULES]):
            raise ValueError('one of the molecules %s not in adsorbent list' %(' '.join(MOLECULES)))
        self.mol2atoms = {row['MOLECULE']: row['ATOMS'].split() for index,row in df_mol.iterrows()}

        if not all([M in list(df_mol['MOLECULE']) for M in MOLECULES]):
            raise ValueError("The molecules mentioned are not supported by the code, see the data directory in the source directory")
        self.MOLECULES = MOLECULES

        df_structures = pd.read_csv(structures_file, encoding='utf-8')
        df_structures = df_structures[['Structures']]
        df_structures['STRUCTURE_NAME'] = df_structures['Structures'].str.replace('.cif','')
        
        if type_ in SIMULATION_TYPES['INFO']:
            self.data = df_structures[['STRUCTURE_NAME','Structures']].to_records(index=False)
        elif type_ in SIMULATION_TYPES['RASPA2']+SIMULATION_TYPES['ZEO++']:
            df_info = pd.read_csv(os.path.join(SOURCE_DIR, "../data/info.csv"), encoding='utf-8') 
            df = pd.merge(df_structures[['STRUCTURE_NAME']], df_info[['STRUCTURE_NAME','UnitCell','Volume [nm^3]']],how="inner", on="STRUCTURE_NAME")
            df = df[df['Volume [nm^3]'] <= Threshold_volume]    
            if type_ in SIMULATION_TYPES['RASPA2']:
                self.data = df[['STRUCTURE_NAME','UnitCell']].to_records(index=False)
            else:
                df['ProbeRadius'] = probe_radius
                self.data = df[['STRUCTURE_NAME','ProbeRadius']].to_records(index=False)

        if setup==True:
            print_every = N_cycles // 10
            init_cycles = min(cycles // 2, 10000)
            self.generate_files(OUTPUT_PATH, type_, FORCE_FIELD=force_field, N_cycles=cycles, N_print=print_every, N_init=init_cycles)


    def generate_files(self, path_to_work, type_, **kwargs):
        """Generate the files need for the simulations

        Args:
            path_to_work (str): path to the working directory, where the simulations occur

        """

        # create input
        path_to_input = os.path.join(SOURCE_DIR, "../Raspa_screening_templates/INPUT_%s"%type_)
        generate(path_to_input, **kwargs)


        if type_ in SIMULATION_TYPES['RASPA2']:
            path_to_Scripts = os.path.join(path_to_work, 'Scripts')
            if not os.path.exists(path_to_Scripts):
                os.mkdir(path_to_Scripts)
            copy_file(SOURCE_DIR,"../Raspa_screening_templates/INPUT_%s"%type_,path_to_work,"INPUT")
            copy_file(SOURCE_DIR,"../Raspa_screening_templates/run",path_to_work,"run")
            if type_ != 'grid':
                copy_file(SOURCE_DIR,"../Raspa_screening_templates/data_%s.sh"%type_,path_to_work,"data.sh")
        elif type_ == "info":
            path_to_Scripts = os.path.join(path_to_work, 'Scripts')
            if not os.path.exists(path_to_Scripts):
                os.mkdir(path_to_Scripts)
            copy_file(SOURCE_DIR,"../Raspa_screening_templates/INPUT_%s"%type_,path_to_work,"INPUT")
            copy_file(SOURCE_DIR,"../Raspa_screening_templates/run_%s"%type_,path_to_work,"run")
        else: # zeo++ modules
            path_to_output = os.path.join(path_to_work, 'Output')
            if not os.path.exists(path_to_output):
                os.mkdir(path_to_output)
            copy_file(SOURCE_DIR,"../Zeo++_screening_templates/run_%s"%type_,path_to_work,"run")

        self.path_to_run = os.path.join(path_to_work,"run")
        path_to_data = os.path.join(path_to_work,"data.sh")
        path_to_INPUT = os.path.join(path_to_work, 'INPUT')
        sed("FORCE_FIELD",FORCE_FIELD,path_to_INPUT)
        if CYCLES:
            sed("N_cycles",CYCLES,path_to_INPUT)
            sed("NCYCLES",CYCLES,self.path_to_run)

        ATOMS = []
        for molecule in MOLECULES:
            atoms = self.mol2atoms[molecule]
            ATOMS += atoms

        sed("N_ATOMS", len(ATOMS), path_to_INPUT)
        sed("ATOMS",' '.join(ATOMS),path_to_INPUT)
        sed("PRESSURES",' '.join(PRESSURES),path_to_INPUT)
        if self.type_ == 'coad':
            sed("MOLECULE",MOLECULES[0],path_to_INPUT)
            sed("MOLECULE",MOLECULES[0],path_to_data)
        elif len(MOLECULES) == 2:
            Y = [round(float(c)/100,6) for c in COMPOSITION]
            sed("MOLECULE_1",MOLECULES[0],path_to_INPUT)
            sed("Y_1",Y[0],path_to_INPUT)
            sed("MOLECULE_2",MOLECULES[1],path_to_INPUT)
            sed("Y_2",Y[1],path_to_INPUT)


    @staticmethod
    def generate(path, **replace_string):
        """Read an input file from Raspa_screening_templates and replace some key words

        Args:
            path (str): path to the INPUT or run or data file (Raspa_screening_templates)

        """
        if not os.path.exists(path):
            return None
        else:
            generated_file = open(path,"r")
            for key, value in replace_string.items():
                generated_file.replace(key, value)
            return generated_file


    @staticmethod
    def write_file(generated_file, outfile_path):
        """Write 

        Args:
            generated_file (str): content of the generated file
            outfile_path (str): path to the output file 

        """

        with open(outfile_path, "w") as outfile:
            outfile.write(generated_file)
            outfile_path.close()


    def run(self, inputs):
        """A module to run one simulation on a given node according to the current process id"""

        FRAMEWORK_NAME,UNITCELL = inputs
        if len(self.NODES) == 0:
            command = "bash %s %s \"%s\""%(self.path_to_run,FRAMEWORK_NAME,UNITCELL)
        else:
            worker = int(mp.current_process()._identity[0])
            nnode = len(self.NODES)
            index = (worker-1)%nnode
            command = "ssh %s \"bash %s %s \\\"%s\\\"\""%(self.NODES[index],self.path_to_run,FRAMEWORK_NAME,UNITCELL)
        return os.system(command)

    def mp_run(self):
        """Using multiprocessing, runs in parallel the simulations"""

        t0 = time()
        print("%s simulation of %s"%(self.type_,' '.join(self.MOLECULES)))
        with mp.Pool(processes=self.nprocs) as p:
            p.map(self.run, [(FRAMEWORK_NAME,UNITCELL) for FRAMEWORK_NAME,UNITCELL in self.data])
        print("SIMULATIONS COMPLETED after %.2f hour(s)"%((time()-t0)/3600))
