import sys
import os
import csv
from time import time
from textwrap import dedent
import multiprocessing as mp

import numpy as np
import pandas as pd

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SOURCE_DIR)
MATSCREEN = os.path.dirname(SOURCE_DIR)

# TODO
# Debug the file writing of Home type simulations (bug when several processes write on a same file)

class Screening():
    def __init__(self, structures_file, procs_per_node, nprocs, type_='grid', force_field="UFF", MOLECULES=['xenon','krypton'], composition=None,
    pressures=[101300], temperature=298.0, cycles=2000, cutoff=12.0, rejection=0.85, probe_radius=1.2, Threshold_volume=20, OUTPUT_PATH=".", RESTART="no"):
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
            glost_list        (bool): Print a list of commands to be run by glost_launch binary
                                      (see https://github.com/cea-hpc/glost.git)

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
        self.SIMULATION_TYPES = {"RASPA2" : ['grid', 'ads', 'coad', 'ent', 'widom', 'widom_nogrid', 'vf', 'sp', 'diffusion','sa'],
                    "INFO"   : ['info'],
                    "ZEO++"  : ['surface', 'volume', 'pore', 'channel', 'voronoi', 'block'],
                    "HOME"   : ['sample', 'surface_sample', 'findsym'],
                    "CPP"    : ["raess", "csurface", "csurface_spiral", "csurface_radius", "csurface_acc","csurface_sa"]
                   }
        try:
            self.NODES = os.environ['NODES']
        except:
            self.NODES = ''
        self.NODES = self.NODES.split()
        if not os.path.exists(os.path.join(MATSCREEN, "set_environment")):
            raise FileNotFoundError("Please check your set_environement file in %s"%MATSCREEN)
        if len(self.NODES)==0:
            if procs_per_node > nprocs:
                raise ValueError('More processes at a time than processors available!!!')
        else:
            if len(self.NODES)*procs_per_node > nprocs:
                raise ValueError('More processes at a time than processors available!!!')
        self.nprocs = nprocs
        self.type_ = type_
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

        molecules_path = os.path.join(MATSCREEN, 'data/molecules.csv')
        df_mol = pd.read_csv(molecules_path, encoding='utf-8')
        if not all([molecule in list(df_mol['MOLECULE']) for molecule in MOLECULES]):
            raise ValueError('One of the molecules %s not in adsorbent list defined in %s'%(' '.join(MOLECULES), molecules_path))
        mol2atoms = {row['MOLECULE']: row['ATOMS'] for index,row in df_mol.iterrows()}

        print_every = cycles // 10
        init_cycles = min(cycles // 2, 10000)
        ATOMS = ' '.join([mol2atoms[molecule] for molecule in MOLECULES])
        N_ATOMS = len(ATOMS.split())
        molecule_dict = {}
        for i in range(len(mole_fraction)):
            molecule_dict[MOLECULES[i]] = mole_fraction[i]

        current_directory = os.getcwd()
        self.OUTPUT_PATH = os.path.join(current_directory, OUTPUT_PATH)
        self.n_sample = cycles
        self.acc_coeff = probe_radius
        self.rej_coeff = rejection

        self.generate_files(self.OUTPUT_PATH, type_, molecule_dict=molecule_dict,
                                                     FORCE_FIELD=force_field,
                                                     N_cycles=cycles,
                                                     N_print=print_every,
                                                     N_init=init_cycles,
                                                     CUTOFF=cutoff,
                                                     PRESSURES=' '.join(pressures),
                                                     TEMPERATURE=temperature,
                                                     N_ATOMS=N_ATOMS,
                                                     ATOMS=ATOMS,
                                                     MOLECULE=MOLECULES[0],
                                                     RESTART=RESTART,
                                                     TIMESTEP=probe_radius,
                                                     REJECT=rejection,
                                                     PATH=self.OUTPUT_PATH,
                                                     )

        df_structures = pd.read_csv(os.path.join(current_directory, structures_file), encoding='utf-8')
        df_structures = df_structures[['Structures']]
        df_structures['STRUCTURE_NAME'] = df_structures['Structures'].str.replace('.cif','', regex=False)

        self.home = False
        self.forcefield = force_field
        if type_ in self.SIMULATION_TYPES['INFO']:
            self.data = df_structures[['STRUCTURE_NAME','Structures']].to_records(index=False)
        elif type_ in self.SIMULATION_TYPES['RASPA2']+self.SIMULATION_TYPES['ZEO++']+self.SIMULATION_TYPES['HOME']+self.SIMULATION_TYPES['CPP']:
            df_info = pd.read_csv(os.path.join(MATSCREEN, "data/info.csv"), encoding='utf-8').drop_duplicates(subset=['STRUCTURE_NAME'])
            df = pd.merge(df_structures[['STRUCTURE_NAME']], df_info[['STRUCTURE_NAME', 'UnitCell', 'Volume [nm^3]','unit vector a', 'unit vector b', 'unit vector c']],how="left", on="STRUCTURE_NAME")
            df['UnitCell'] = df['UnitCell'].fillna("1 1 1")
            df['Volume [nm^3]'] = df['Volume [nm^3]'].fillna(Threshold_volume)
            df = df[df['Volume [nm^3]'] <= Threshold_volume]
            if type_ in self.SIMULATION_TYPES['RASPA2']+self.SIMULATION_TYPES['CPP']:
                self.data = df[['STRUCTURE_NAME','UnitCell']].to_records(index=False)
            elif type_ in self.SIMULATION_TYPES['ZEO++']:
                df['ProbeRadius'] = probe_radius
                self.data = df[['STRUCTURE_NAME','ProbeRadius']].to_records(index=False)
            elif type_ in self.SIMULATION_TYPES['HOME']:
                self.atoms = '|'.join([mol2atoms[molecule] for molecule in MOLECULES])
                self.temperature = temperature
                self.cutoff = float(cutoff)
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
            path_to_INPUT = os.path.join(MATSCREEN, "Raspa_screening_templates/INPUT_%s"%type_)
            INPUT_file = self.generate(path_to_INPUT, **kwargs)
            if type_ == 'coad':
                index = 0
                for key, value in molecule_dict.items():
                    INPUT_file += dedent("""
                    Component %d MoleculeName                     %s
                                MoleculeDefinition               TraPPE
                                TranslationProbability           0.5
                                IdentityChangeProbability        1.0
                                NumberOfIdentityChanges          2
                                IdentityChangesList              0 1
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
                RUN_file = open(os.path.join(MATSCREEN, "Raspa_screening_templates/run_%s"%type_), "r").read()
            if type_ == "diffusion":
                RUN_file = self.generate(os.path.join(MATSCREEN, "Raspa_screening_templates/run_diffusion"), **kwargs)
                os.system("mkdir %s/Output"%path_to_work)
            else:
                RUN_file = open(os.path.join(MATSCREEN, "Raspa_screening_templates/run"), "r").read()

            self.path_to_run = os.path.join(path_to_work, "run")
            self.write_file(RUN_file, self.path_to_run)
            if type_ != 'grid':
                DATA_file = self.generate(os.path.join(MATSCREEN, "Raspa_screening_templates/data_%s.sh"%type_), **kwargs)
                self.write_file(DATA_file, os.path.join(path_to_work,"data.sh"))
                if type_ == 'info':
                    merge_file = self.generate(os.path.join(MATSCREEN, "Raspa_screening_templates/merge_info.py"), **kwargs)
                    self.write_file(merge_file, os.path.join(path_to_work,"merge_info.py"))

        elif type_ in self.SIMULATION_TYPES["ZEO++"]:
            path_to_Output = os.path.join(path_to_work, 'Output')
            if not os.path.exists(path_to_Output):
                os.mkdir(path_to_Output)
            if type_ == "voronoi":
                if not os.path.exists(os.path.join(path_to_work, 'Coordinates')):
                    os.mkdir(os.path.join(path_to_work, 'Coordinates'))
                os.system("cp %s %s"%(os.path.join(MATSCREEN, "Zeo++_screening_templates/extract_vertex.py"),path_to_work))
            RUN_file = self.generate(os.path.join(MATSCREEN, "Zeo++_screening_templates/run_%s"%type_), **kwargs)
            self.path_to_run = os.path.join(path_to_work,"run")
            self.write_file(RUN_file, self.path_to_run)

        elif type_ in self.SIMULATION_TYPES["HOME"]:
            self.path_to_run = os.path.join(path_to_work,"run.py")
            os.system("cp %s %s"%(os.path.join(MATSCREEN, "Home_screening_templates/run_%s.py"%type_),self.path_to_run))
            if type_ == "findsym":
                path_to_Output = os.path.join(path_to_work, 'Output')
                if not os.path.exists(path_to_Output):
                    os.mkdir(path_to_Output)
            else:
                pd.DataFrame(columns={"Structure_name":[], "Adsorbent_name":[], "Acessible_average_energy":[], "Minimum_energy":[], "Boltzmann_average_energy":[], "Henry_coeff":[]}).to_csv('home_output.csv',index=False)
                open(os.path.join(path_to_work, '.output_written.tmp'), 'a')

        elif type_ in self.SIMULATION_TYPES["CPP"]:
            self.path_to_run = os.path.join(path_to_work,"run.sh")
            RUN_file = self.generate(os.path.join(MATSCREEN, "Cpp_screening_templates/run_%s.sh"%type_), **kwargs)
            self.write_file(RUN_file, self.path_to_run)
            if type_ in ['raess']:
                pd.DataFrame(columns={"Structure_name":[], "Enthalpy_surface_kjmol":[], "Henry_coeff_molkgPa":[], "ASA_m2_cm3":[], "time":[]}).to_csv('cpp_output_%s_%s_%s.csv'%(self.n_sample,self.rej_coeff,self.acc_coeff),index=False)
            elif type_ in ['csurface_acc','csurface_sa']:
                pd.DataFrame(columns={"Structure_name":[], "Enthalpy_surface_kjmol":[], "Henry_coeff_molkgPa":[], "time":[]}).to_csv('cpp_output_%s.csv'%(self.acc_coeff),index=False)
            else:
                pd.DataFrame(columns={"Structure_name":[], "Enthalpy_surface_kjmol":[], "Henry_coeff_molkgPa":[], "time":[]}).to_csv('cpp_output_%s.csv'%(self.n_sample),index=False)

        #os.system("bash %s %s"%(os.path.join(MATSCREEN, "copy_env.sh"), os.path.join(path_to_work, "set_environment")))


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
            command = "bash %s %s \"%s\" %s"%(self.path_to_run,FRAMEWORK_NAME,UNITCELL,self.forcefield)
        else:
            worker = int(mp.current_process()._identity[0])
            nnode = len(self.NODES)
            index = (worker-1)%nnode
            HOST = self.NODES[index]
            command = "ssh %s \"bash %s %s \\\"%s\\\" %s\""%(HOST,self.path_to_run,FRAMEWORK_NAME,UNITCELL,self.forcefield)
        os.system(command)
        output_dict = {'Structures':[FRAMEWORK_NAME], "CPU_time (s)":[int(time()-t0)]}
        pd.DataFrame(output_dict).to_csv(os.path.join(self.OUTPUT_PATH,"time.csv"),mode="a",index=False,header=False)


    def run_home(self, inputs):
        """Runs calculations using a homemade python algorithm based on pymatgen
        """

        t0 = time()
        structure_name, supercell_wrap = inputs
        if self.type_ == "surface_sample":
            supercell_wrap = self.n_sample
        if len(self.NODES) == 0:
            command = "%s %s %s %s %s %s \"%s\" \"%s\" "%(sys.executable, self.path_to_run, structure_name, self.cutoff, self.forcefield, self.temperature, self.atoms, supercell_wrap)
            print(command)
        else:
            worker = int(mp.current_process()._identity[0])
            nnode = len(self.NODES)
            index = (worker-1)%nnode
            HOST = self.NODES[index]
            command = "ssh %s \"%s %s \\\"%s\\\" %s %s %s %s \\\"%s\\\" \""%(HOST, sys.executable, self.path_to_run, self.atoms, self.forcefield, self.temperature, self.cutoff, structure_name, supercell_wrap)
        os.system(command)
        output_dict = {'Structures':[structure_name], "CPU_time (s)":[int(time()-t0)]}
        pd.DataFrame(output_dict).to_csv(os.path.join(self.OUTPUT_PATH,"time.csv"),mode="a",index=False,header=False)

    def glost_list(self):
        """ Print out glost list for mpirun"""
        df = pd.DataFrame.from_records(self.data)
        struc = df.iloc[:,0]
        opt = df.iloc[:,1]
        if self.home:
            command = "%s %s \"%s\" %s %s %s "%(sys.executable, self.path_to_run, self.atoms, self.forcefield, self.temperature, self.cutoff) + struc + " \"" + opt + "\""
        else:
            command = "bash %s "%(self.path_to_run) + struc + " \"" + opt + "\""
        command.to_frame().to_csv(os.path.join(self.OUTPUT_PATH,"glost.list"),index=False,header=False, quoting=csv.QUOTE_NONE, quotechar='')


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
