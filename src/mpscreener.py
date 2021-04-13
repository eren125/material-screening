import sys
import os
from time import time
from textwrap import dedent
import multiprocessing as mp

import numpy as np
import pandas as pd

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SOURCE_DIR)

# TODO 


class Screening():
    def __init__(self, structures_file, procs_per_node, nprocs, type_='grid', force_field="UFF", MOLECULES=['xenon','krypton'], composition=None, 
    pressures=[101300], temperature=298.0, cycles=2000, cutoff=12.0, probe_radius=1.2, Threshold_volume=20, OUTPUT_PATH="."):
        """A class for screening purposes using Raspa2 for molecular simulations
        Initialise important variables like the name of the structures to screen and the unitcell associated
        Catch Obvious value errors, incompatible mix of varibles, etc.
        If setup is True, creates file for the simulation INPUT and run mainly (data for post-processing and RestartInitial for single point calculations)

        Args:
            structures_file    (str): relative path to the csv file containing all the structures, 
                                      the headers must contain "Structures"
            procs_per_node     (int): number of processors available in each node for a given calculator
            nprocs             (int): number of processes to be run at the same time
            type_              (str): type of simulation the user wants to carry out. Currently available:
                                      grid calculation via 'grid', GCMC via 'ads' 'coad', NVT MC via 'ent',
                                      Widom's insertion via 'widom', helium void fraction via 'vf', 
                                      Zeo++ via 'surface' 'volume' 'pore' 'channel' and global information via 'info'
            force_field        (str): force field used for the molecular simulations, 
                                      it must be defined in the Raspa directory 
            MOLECULES         (list): list of molecules (str) to be adsorbed on the materials
            composition       (list): list of mole fractions for the molecules defined in MOLECULES
                                      must have the same length as MOLECULES, and sum must equal to 1
            pressures         (list): list of pressures (float) in pascal to be simulated
            TEMP             (float): temperature in kelvin for the Raspa2 simulations (default = 298.0 K)
            cycles             (int): number of prodution cycles used for the Raspa2 simulations
                                      equilibration cycles are fixed to 10k for now (to improve)                 
            cutoff           (float): sets the van der Waals cutoff in Raspa2 simulation
            probe_radius     (float): radius of the probe in angström considered in Zeo++ simulations
            Threshold_volume (float): threshold at which the volume is considered too big for grid calculations
                                      it consumes too much RAM Memory for the machine currently used
            OUTPUT_PATH        (str): relative path to where the outputs of the simulations will be printed

        self variables:
            NODES              (str): names of the nodes allocated by the supercalculator
                                      if you are running on a local computer 'NODES' does not exist
                                      it will then be set to an empty string by default

        self methods:
            generate_files: generate the input and run files needed for the simulations according to the type_
                            and replace some keywords with the variables given to the function
            generate      : generate an input string according to the file specified
            write_file    : print out a string in a file at the specified output path
            run           : function that takes the framework's name and the smallest unitcell for a 12 angstrom cut-off
                            and runs the corresponding simulation according to the run bash file 
            run_mp        : function that calls `run` multiple times in parallel. mp.Pool distributes the jobs 
                            to $nprocs workers so that every job  
        """

        ### Initialisation of class objects and error catching ###
        self.SIMULATION_TYPES = {"RASPA2" : ['grid', 'ads', 'coad', 'ent', 'widom', 'vf', 'sp'],
                    "INFO"   : ['info'],
                    "ZEO++"  : ['surface', 'volume', 'pore', 'channel', 'voronoi'],
                    "HOME"   : ['sample'] 
                   }
        try:
            self.NODES = os.environ['NODES']
        except:
            self.NODES = ''
        self.NODES = self.NODES.split()
        try:
            os.system("source %s/../set_environment"%SOURCE_DIR)
        except:
            raise FileNotFoundError("Please check your set_environement file in %s"%SOURCE_DIR)
        if len(self.NODES)==0:
            if procs_per_node > nprocs:
                raise ValueError('More processes at a time than proccessors available!!!')
        else:
            if len(self.NODES)*procs_per_node > nprocs:
                raise ValueError('More processes at a time than proccessors available!!!')
        self.nprocs = nprocs

        available_types = sum(self.SIMULATION_TYPES.values(), [])
        if type_ not in available_types:
            raise ValueError(('%s not an option yet. Please choose between: ' +
                ', '.join(['%s']*len(available_types))) % tuple([type_]+available_types))

        if (type_=='coad') and (not composition):
            raise ValueError("Please specify a composition, for example: 20 80")

        if composition:
            mole_fraction = [round(float(c)/100, 4) for c in composition]
            if sum([float(c) for c in composition]) != 100.0:
                raise ValueError("Composition does not add up to 100%%")
            if len(composition) > len(MOLECULES):
                raise ValueError("More composition values than molecules specified")
        else:
            mole_fraction = []

        df_mol = pd.read_csv(os.path.join(SOURCE_DIR, '../data/molecules.csv'), encoding='utf-8')
        if not all([molecule in list(df_mol['MOLECULE']) for molecule in MOLECULES]):
            raise ValueError('one of the molecules %s not in adsorbent list' %(' '.join(MOLECULES)))
        mol2atoms = {row['MOLECULE']: row['ATOMS'] for index,row in df_mol.iterrows()}

        if not all([M in list(df_mol['MOLECULE']) for M in MOLECULES]):
            raise ValueError("The molecules mentioned are not supported by the code, see the data directory in the source directory")

        print_every = cycles // 10
        init_cycles = min(cycles // 2, 10000)
        ATOMS = ' '.join([mol2atoms[molecule] for molecule in MOLECULES])
        N_ATOMS = len(ATOMS.split())
        molecule_dict = {}
        for i in range(len(mole_fraction)):
            molecule_dict[MOLECULES[i]] = mole_fraction[i]
        
        current_directory = os.environ['CURRENTDIR']
        self.OUTPUT_PATH = os.path.join(current_directory, OUTPUT_PATH)

        self.generate_files(self.OUTPUT_PATH, type_, molecule_dict=molecule_dict, FORCE_FIELD=force_field, N_cycles=cycles, N_print=print_every, N_init=init_cycles, 
        CUTOFF=cutoff, PRESSURES=' '.join(pressures), TEMPERATURE=temperature, N_ATOMS=N_ATOMS, ATOMS=ATOMS, MOLECULE=MOLECULES[0])

        df_structures = pd.read_csv(os.path.join(current_directory, structures_file), encoding='utf-8')
        df_structures = df_structures[['Structures']]
        df_structures['STRUCTURE_NAME'] = df_structures['Structures'].str.replace('.cif','', regex=False)
        
        self.home = False
        if type_ in self.SIMULATION_TYPES['INFO']:
            self.data = df_structures[['STRUCTURE_NAME','Structures']].to_records(index=False)
        elif type_ in self.SIMULATION_TYPES['RASPA2']+self.SIMULATION_TYPES['ZEO++']+self.SIMULATION_TYPES['HOME']:
            df_info = pd.read_csv(os.path.join(SOURCE_DIR, "../data/info.csv"), encoding='utf-8') 
            df = pd.merge(df_structures[['STRUCTURE_NAME']], df_info[['STRUCTURE_NAME', 'UnitCell', 'Volume [nm^3]','unit vector a', 'unit vector b', 'unit vector c']],how="inner", on="STRUCTURE_NAME")
            df = df[df['Volume [nm^3]'] <= Threshold_volume]    
            if type_ in self.SIMULATION_TYPES['RASPA2']:
                self.data = df[['STRUCTURE_NAME','UnitCell']].to_records(index=False)
            elif type_ in self.SIMULATION_TYPES['ZEO++']:
                df['ProbeRadius'] = probe_radius
                self.data = df[['STRUCTURE_NAME','ProbeRadius']].to_records(index=False)
            elif type_ in self.SIMULATION_TYPES['HOME']:
                self.atoms = '|'.join([mol2atoms[molecule] for molecule in MOLECULES])
                self.forcefield = force_field
                self.temperature = temperature
                self.cutoff = cutoff
                df['supercell_wrap'] = df['UnitCell'].apply(lambda x: x.replace(' ','|'))
                self.data = df[['STRUCTURE_NAME','supercell_wrap']].to_records(index=False)
                self.home = True

        pd.DataFrame({'Structures':[],"CPU_time (s)":[]}).to_csv(os.path.join(self.OUTPUT_PATH,"time.csv"), index=False)
        print("%s simulation of %s"%(type_,' '.join(MOLECULES)))


    def generate_files(self, path_to_work, type_, molecule_dict={}, **kwargs):
        """Generate the files need for the simulations

        Args:
            path_to_work (str): path to the working directory, where the simulations occur
        """

        if type_ in self.SIMULATION_TYPES['RASPA2']+self.SIMULATION_TYPES["INFO"]:
            path_to_Scripts = os.path.join(path_to_work, 'Scripts')
            if not os.path.exists(path_to_Scripts):
                os.mkdir(path_to_Scripts)
            path_to_INPUT = os.path.join(SOURCE_DIR, "../Raspa_screening_templates/INPUT_%s"%type_)
            INPUT_file = self.generate(path_to_INPUT, **kwargs)
            if type_ == 'coad':
                index = 0
                for key, value in molecule_dict.items():
                    INPUT_file += dedent("""
                    Component %d MoleculeName                     %s
                                MoleculeDefinition               TraPPE
                                TranslationProbability           0.5
                                IdentityChangeProbability        1.0
                                  NumberOfIdentityChanges        2
                                  IdentityChangesList            0 1
                                SwapProbability                  1.0
                                CreateNumberOfMolecules          0
                                MolFraction                      %s     
                    """%(index,key,value))
                    index += 1
            self.write_file(INPUT_file, os.path.join(path_to_work, "INPUT"))
            if type_ == 'sp' or type_ == 'sample':
                if not os.path.exists(os.path.join(path_to_work, "Coordinates")):
                    raise FileNotFoundError("Check that Coordinates/ directory exists and it is loaded with results from Voronoi simulations")
                if not os.path.exists(os.path.join(path_to_work, "RestartInitial/System_0")):
                    os.mkdir(os.path.join(path_to_work, "RestartInitial"))
                    os.mkdir(os.path.join(path_to_work, "RestartInitial/System_0"))
                RUN_file = open(os.path.join(SOURCE_DIR, "../Raspa_screening_templates/run_%s"%type_), "r").read()
            else:
                RUN_file = open(os.path.join(SOURCE_DIR, "../Raspa_screening_templates/run"), "r").read()
            self.path_to_run = os.path.join(path_to_work,"run")
            self.write_file(RUN_file, self.path_to_run)
            if type_ != 'grid':
                DATA_file = open(os.path.join(SOURCE_DIR,"../Raspa_screening_templates/data_%s.sh"%type_), "r").read()
                self.write_file(DATA_file, os.path.join(path_to_work,"data.sh"))
                if type_ == 'info':
                    os.system("cp %s %s"%(os.path.join(SOURCE_DIR, "../Raspa_screening_templates/merge_info.py"), os.path.join(path_to_work,"merge_info.py")))

        elif type_ in self.SIMULATION_TYPES["ZEO++"]:
            path_to_Output = os.path.join(path_to_work, 'Output')
            if not os.path.exists(path_to_Output):
                os.mkdir(path_to_Output)
            if type_ == "voronoi":
                if not os.path.exists(os.path.join(path_to_work, 'Coordinates')):
                    os.mkdir(os.path.join(path_to_work, 'Coordinates'))
                os.system("cp %s %s"%(os.path.join(SOURCE_DIR,"../Zeo++_screening_templates/extract_vertex.py"),path_to_work))
            RUN_file = open(os.path.join(SOURCE_DIR,"../Zeo++_screening_templates/run_%s"%type_), "r").read()
            self.path_to_run = os.path.join(path_to_work,"run")
            self.write_file(RUN_file, self.path_to_run)

        elif type_ in self.SIMULATION_TYPES["HOME"]:
            os.system("cp %s %s"%(os.path.join(SOURCE_DIR,"../Home_screening_templates/run.py"),path_to_work))
            self.path_to_run = os.path.join(path_to_work,"run.py")
            pd.DataFrame(columns={"Structure_name":[], "Adsorbent_name":[], "Acessible_average_energy":[], "Minimum_energy":[], "Boltzmann_average_energy":[]}).to_csv('home_output.csv',index=False)


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
        """Write a given string in the path specified

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

        t0 = time()
        FRAMEWORK_NAME,UNITCELL = inputs
        if len(self.NODES) == 0:
            command = "bash %s %s \"%s\""%(self.path_to_run,FRAMEWORK_NAME,UNITCELL)
        else:
            worker = int(mp.current_process()._identity[0])
            nnode = len(self.NODES)
            index = (worker-1)%nnode
            HOST = self.NODES[index]
            command = "ssh %s \"bash %s %s \\\"%s\\\"\""%(HOST,self.path_to_run,FRAMEWORK_NAME,UNITCELL)
        os.system(command)
        output_dict = {'Structures':[FRAMEWORK_NAME], "CPU_time (s)":[int(time()-t0)]}
        pd.DataFrame(output_dict).to_csv(os.path.join(self.OUTPUT_PATH,"time.csv"),mode="a",index=False,header=False)


    def run_home(self, inputs):
        """Runs calculations using a homemade python algorithm based on pymatgen
        """

        t0 = time()
        structure_name, supercell_wrap = inputs
        if len(self.NODES) == 0:
            command = "python3 %s \"%s\" %s %s %s %s \"%s\" "%(self.path_to_run, self.atoms, self.forcefield, self.temperature, self.cutoff, structure_name, supercell_wrap)
            print(command)
        else:
            worker = int(mp.current_process()._identity[0])
            nnode = len(self.NODES)
            index = (worker-1)%nnode
            HOST = self.NODES[index]
            command = "ssh %s \"python3 %s \\\"%s\\\" %s %s %s %s \\\"%s\\\" \""%(HOST,self.path_to_run, self.atoms, self.forcefield, self.temperature, self.cutoff, structure_name, supercell_wrap)
        os.system(command)
        output_dict = {'Structures':[structure_name], "CPU_time (s)":[int(time()-t0)]}
        pd.DataFrame(output_dict).to_csv(os.path.join(self.OUTPUT_PATH,"time.csv"),mode="a",index=False,header=False)


    def mp_run(self):
        """Using multiprocessing, runs in parallel the simulations
        """
        t0 = time()

        if self.home:
            with mp.Pool(processes=self.nprocs) as p:
                p.map(self.run_home, [(structure_name, supercell_wrap) for structure_name, supercell_wrap in self.data])
        else: # homemade simulations
            with mp.Pool(processes=self.nprocs) as p:
                p.map(self.run, [(FRAMEWORK_NAME,UNITCELL) for FRAMEWORK_NAME,UNITCELL in self.data])

        print("SIMULATIONS COMPLETED after %.2f hour(s)"%((time()-t0)/3600))
