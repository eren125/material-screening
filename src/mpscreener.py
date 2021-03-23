import sys
import os
SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SOURCE_DIR)

from time import time
import numpy as np
import pandas as pd

import multiprocessing as mp

SIMULATION_TYPES = {"RASPA2" : ['grid', 'ads', 'coad', 'ent', 'widom', 'vf', 'point'],
                    "INFO"   : ['info'],
                    "ZEO++"  : ['surface', 'volume', 'pore', 'channel', 'voronoi']
                   }


# TODO
# Finish coad
# Add timer inside run files and echo in a file to create a database of time spent per sim // 
# can get it from output but more accurate (including preprocessing time)
# Add single point simulation logic
# Chain several sp sim
# Move procs and procs_per_node to mp_run?
# ADD cutoff as a variable for raspa simulations

class Screening():
    def __init__(self, structures_file, force_field, MOLECULES, nprocs, pressures=[101300], temperature=298.0, cutoff=12, probe_radius=1.2, 
    Threshold_volume=20, procs_per_node=48, type_='grid', OUTPUT_PATH=".", composition=None, positions=None, cycles=2000, setup=True):
        """A class for screening purposes using Raspa2 for molecular simulations
        (the path of the Raspa directory have to be exported as en environement variable)

        Args:
            structures_file    (str): relative path to the csv file containing all the structures, 
                                      the headers must contain "Structures"
            force_field        (str): force field used for the molecular simulations, 
                                      it must be defined in the Raspa directory 
            MOLECULES         (list): list of molecules (str) to be adsorbed on the materials
            pressures         (list): list of pressures (float) in pascal to be simulated
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
            composition       (list): list of mole fractions for the molecules defined in MOLECULES
                                      must have the same length as MOLECULES, and sum must equal to 1
            cycles             (int): number of prodution cycles used for the Raspa2 simulations
                                      equilibration cycles are fixed to 10k for now (to improve)
            setup             (bool): boolean variable to determine whether to copy the input files 
                                      (obsolete in next version: the input will be parsed)

        self variables:
            NODES              (str): names of the nodes allocated by the supercalculator
                                      if you are running on a local computer 'NODES' does not exist
                                      it will then be set to an empty string by default

        self methods:
            generate_files
            generate
            write_file
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

        available_types = []
        for key in SIMULATION_TYPES.keys():
            available_types += SIMULATION_TYPES[key]
        if type_ not in available_types:
            raise ValueError(('%s not an option yet. Please choose between: ' +
                ', '.join(['%s']*len(available_types))) % tuple([type_]+available_types))

        if (type_=='point') and (not positions):
            raise ValueError("Please specify the path to a csv file with positions. See examples in data/positions_sample.csv.")

        if (type_=='coad') and (not composition):
            raise ValueError("Please specify a composition, for example: 20 80")

        if composition:
            if sum([float(c) for c in composition]) != 100.0:
                raise ValueError("Composition does not add up to 100%%")

        df_mol = pd.read_csv(os.path.join(SOURCE_DIR, '../data/molecules.csv'), encoding='utf-8')
        if not all([molecule in list(df_mol['MOLECULE']) for molecule in MOLECULES]):
            raise ValueError('one of the molecules %s not in adsorbent list' %(' '.join(MOLECULES)))
        mol2atoms = {row['MOLECULE']: row['ATOMS'] for index,row in df_mol.iterrows()}

        if not all([M in list(df_mol['MOLECULE']) for M in MOLECULES]):
            raise ValueError("The molecules mentioned are not supported by the code, see the data directory in the source directory")

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
            print_every = cycles // 10
            init_cycles = min(cycles // 2, 10000)
            ATOMS = ' '.join([mol2atoms[molecule] for molecule in MOLECULES])
            N_ATOMS = len(ATOMS.split())

            self.generate_files(OUTPUT_PATH, type_, FORCE_FIELD=force_field, N_cycles=cycles, N_print=print_every, N_init=init_cycles, 
            CUTOFF=cutoff, PRESSURES=' '.join(pressures), TEMPERATURE=temperature, N_ATOMS=N_ATOMS, ATOMS=ATOMS, MOLECULE=MOLECULES[0])

            if type_ == "coad":
                # loop over mole fractions and molecules
                s = """
                Component 0 MoleculeName                     MOLECULE_1
                            MoleculeDefinition               TraPPE
                            TranslationProbability           0.5
                            IdentityChangeProbability        1.0
                              NumberOfIdentityChanges        2
                              IdentityChangesList            0 1
                            SwapProbability                  1.0
                            CreateNumberOfMolecules          0
                            MolFraction                      Y_1


                Component 1 MoleculeName                     MOLECULE_2
                            MoleculeDefinition               TraPPE
                            TranslationProbability           0.5
                            IdentityChangeProbability        1.0
                              NumberOfIdentityChanges        2
                              IdentityChangesList            0 1
                            SwapProbability                  1.0
                            CreateNumberOfMolecules          0
                            MolFraction                      Y_2
                """

        print("%s simulation of %s"%(type_,' '.join(MOLECULES)))

    def generate_files(self, path_to_work, type_, **kwargs):
        """Generate the files need for the simulations

        Args:
            path_to_work (str): path to the working directory, where the simulations occur
        """

        # create input
        if type_ in SIMULATION_TYPES['RASPA2']+SIMULATION_TYPES["INFO"]:
            path_to_Scripts = os.path.join(path_to_work, 'Scripts')
            if not os.path.exists(path_to_Scripts):
                os.mkdir(path_to_Scripts)
            path_to_INPUT = os.path.join(SOURCE_DIR, "../Raspa_screening_templates/INPUT_%s"%type_)
            INPUT_file = self.generate(path_to_INPUT, **kwargs)
            self.write_file(INPUT_file, os.path.join(path_to_work, "INPUT"))
            RUN_file = open(os.path.join(SOURCE_DIR, "../Raspa_screening_templates/run"), "r").read()
            self.path_to_run = os.path.join(path_to_work,"run")
            self.write_file(RUN_file, self.path_to_run)
            if type_ != 'grid':
                DATA_file = open(os.path.join(SOURCE_DIR,"../Raspa_screening_templates/data_%s.sh"%type_), "r").read()
                self.write_file(DATA_file, os.path.join(path_to_work,"data.sh"))
        elif type_ in SIMULATION_TYPES["ZEO++"]:
            path_to_Output = os.path.join(path_to_work, 'Output')
            if not os.path.exists(path_to_Output):
                os.mkdir(path_to_Output)
            RUN_file = open(os.path.join(SOURCE_DIR,"../Zeo++_screening_templates/run_%s"%type_), "r").read()
            self.path_to_run = os.path.join(path_to_work,"run")
            self.write_file(RUN_file, self.path_to_run)


    @staticmethod
    def generate(path, **replace_string):
        """Read an input file from Raspa_screening_templates and replace some key words

        Args:
            path (str): path to the INPUT or run or data file (Raspa_screening_templates)
        """
        if not os.path.exists(path):
            return None
        else:
            generated_file = open(path,"r").read()
            for key, value in replace_string.items():
                generated_file = generated_file.replace(str(key), str(value))
            print(generated_file)
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
            outfile.close()


    def run(self, inputs):
        """A module to run one simulation on a given node according to the current process id
        """

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
        """Using multiprocessing, runs in parallel the simulations
        """

        t0 = time()
        with mp.Pool(processes=self.nprocs) as p:
            p.map(self.run, [(FRAMEWORK_NAME,UNITCELL) for FRAMEWORK_NAME,UNITCELL in self.data])
        print("SIMULATIONS COMPLETED after %.2f hour(s)"%((time()-t0)/3600))
