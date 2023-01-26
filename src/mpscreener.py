import sys
import os
import csv
from time import time
from textwrap import dedent
import multiprocessing as mp

import numpy as np
import pandas as pd
from pathlib import Path

import shlex

SOURCE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SOURCE_DIR)
MATSCREEN = os.path.dirname(SOURCE_DIR)

# TODO
# Debug the file writing of Home type simulations (bug when several processes write on a same file)

class Screening():
    def __init__(self, structures_file, procs_per_node, nprocs, type_='grid', force_field="UFF",
                 MOLECULES=['xenon','krypton'], composition=None, pressures=[101300], temperatures=[298.0],
                 cycles=2000, EwaldPrecision=1e-6, cutoff=12.0, rejection=0.85, probe_radius=1.2,
                 Threshold_volume=0, OUTPUT_PATH=".", RESTART=False, N_init=-1, print_every=-1, MOVIE=False, EXTRA="",
                 SKIPDONE=False, PressureParallelism='parallel'):
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
                                      grid calculation via 'grid', GCMC via 'ads' or 'coad', parallel tempering via 'pt',
                                      NVT MC via 'ent', Widom's insertion via 'widom', helium void fraction via 'vf',
                                      Zeo++ via 'surface' 'volume' 'pore' 'channel', global information via 'info'.
            force_field        (str): force field used for the molecular simulations,
                                      it must be defined in the Raspa directory
            MOLECULES         (list): list of molecules (str) to be adsorbed on the materials
            composition       (list): list of mole fractions for the molecules defined in MOLECULES
                                      must have the same length as MOLECULES, and sum must equal to 1
            pressures         (list): list of pressures (float) in pascal to be simulated
            temperatures      (list): list of temperatures in Kelvin for the Raspa2 simulations (default = [298.0] K)
            cycles             (int): number of prodution cycles used for the Raspa2 simulations
            N_init             (int): number of initialization cycle (default is min(cycles//2, 10000),
                                      or 10 times less if RESTART is set)
            cutoff           (float): sets the van der Waals cutoff in Raspa2 simulation
            probe_radius     (float): radius of the probe in angström considered in Zeo++ simulations
            EwaldPrecision   (float): precision for Ewald summation, or 0 if no Ewald is used
            Threshold_volume (float): threshold at which the volume is considered too big for grid calculations
                                      it consumes too much RAM Memory for the machine currently used
            OUTPUT_PATH        (str): relative path to where the outputs of the simulations will be printed
            glost_list        (bool): Print a list of commands to be run by glost_launch binary
                                      (see https://github.com/cea-hpc/glost.git)
            RESTART           (bool): specify if the RASPA simulation should be restarted
            SKIPDONE          (bool): specify whether to check for already completed RASPA simulations and skip them
            MOVIE             (bool): specify if a movie should be exported for each molecule
            EXTRA              (str): extra arguments appended to each framework in Raspa2 INPUT file
            PressureParallelism(str): whether Raspa computations with different pressures should be computed in
                                      sequentially with the default Raspa behaviour ('raspa'), by using the last
                                      Restart file as starting point for each next computation ('sequence'),
                                      or in parallel ('parallel', default)

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
        self.SIMULATION_TYPES = {"RASPA2" : ['grid', 'ads', 'coad', 'ent', 'widom', 'widom_nogrid', 'vf', 'sp', 'diffusion','sa', 'pt', 'samplerestart'],
                    "INFO"   : ['info'],
                    "ZEO++"  : ['surface', 'volume', 'pore', 'channel', 'voronoi', 'block', 'psd','strinfo','volpo','oms'],
                    "HOME"   : ['sample', 'surface_sample', 'findsym'],
                    "CPP"    : ["raess", "cgrid", "barrier", "csurface", "csurface_spiral", "csurface_radius", "csurface_acc","csurface_sa"]
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
        if print_every == -1:
            print_every = min(cycles//10,100)

        if N_init == -1:
            N_init = min(cycles//2, 10000)//(9*RESTART+1)
        self.N_init = N_init
        assert PressureParallelism == "parallel" or PressureParallelism == "sequence" or PressureParallelism == "raspa"
        self.PressureParallelism = PressureParallelism
        self.pressures = pressures

        mole_fraction = []
        if type_ == 'coad':
            if not composition:
                raise ValueError("Please specify a composition, for example: 20 80")
            mole_fraction = [round(float(c)/100, 4) for c in composition]
            if sum([float(c) for c in composition]) != 100.0:
                raise ValueError("Composition does not add up to 100%%")
            if len(composition) > len(MOLECULES):
                raise ValueError("More composition values than molecules specified")
        elif type_ == 'ads' or type_ == 'samplerestart':
            if not composition:
                print("WARNING: no composition specified! All components are assumed to be molecules (not cations).\n")
                composition = [0 for _ in MOLECULES]
            elif len(composition) < len(MOLECULES):
                print("WARNING: some components do not have corresponding nature (cation or molecules). They will be considered as molecules.")
                composition.extend(0 for _ in range(len(MOLECULES) - len(composition)))
            mole_fraction = [bool(int(x)) for x in composition]
            if type_ != 'samplerestart' and any(mole_fraction):
                RESTART = True
        elif type_ == 'pt':
            mole_fraction = [True for _ in MOLECULES]

        self.RESTART = RESTART
        self.SKIPDONE = SKIPDONE

        molecules_path = os.path.join(MATSCREEN, 'data/molecules.csv')
        df_mol = pd.read_csv(molecules_path, encoding='utf-8')
        if not all([molecule in list(df_mol['MOLECULE']) for molecule in MOLECULES]):
            raise ValueError('One of the molecules %s not in adsorbent list defined in %s'%(' '.join(MOLECULES), molecules_path))
        mol2atoms = {row['MOLECULE']: row['ATOMS'] for index,row in df_mol.iterrows()}

        ATOMS = ' '.join([mol2atoms[molecule] for molecule in MOLECULES])
        N_ATOMS = len(ATOMS.split())
        molecule_dict = {}
        for i in range(len(mole_fraction)):
            molecule_dict[MOLECULES[i]] = mole_fraction[i]

        current_directory = os.getcwd()
        self.OUTPUT_PATH = os.path.join(current_directory, OUTPUT_PATH)
        self.n_sample = cycles
        self.acc_coeff = probe_radius
        EWALD = "Ewald\nEwaldPrecision "+str(EwaldPrecision) if EwaldPrecision else "None"
        self.movie = MOVIE
        self.temperatures = temperatures
        self.rej_coeff = rejection

        self.generate_files(self.OUTPUT_PATH, type_, mol2atoms, molecule_dict=molecule_dict,
                                                     FORCE_FIELD=force_field,
                                                     N_cycles=cycles,
                                                     N_print=print_every,
                                                     CUTOFF=cutoff,
                                                     EWALD=EWALD,
                                                     N_ATOMS=N_ATOMS,
                                                     ATOMS=ATOMS,
                                                     MOLECULE=MOLECULES[0],
                                                     TIMESTEP=probe_radius,
                                                     REJECT=rejection,
                                                     PATH=self.OUTPUT_PATH,
                                                     EXTRA=EXTRA,
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
            if cutoff != 12 :
                df_info['a'] = ((cutoff*2)/df_info['x_box [A]'].astype(float)).astype(int) + 1
                df_info['b'] = ((cutoff*2)/df_info['y_box [A]'].astype(float)).astype(int) + 1
                df_info['c'] = ((cutoff*2)/df_info['z_box [A]'].astype(float)).astype(int) + 1
                df_info['UnitCell'] = df_info['a'].astype(str) + ' ' + df_info['b'].astype(str) + ' ' + df_info['c'].astype(str)
            df = pd.merge(df_structures[['STRUCTURE_NAME']], df_info[['STRUCTURE_NAME', 'UnitCell', 'Volume [nm^3]']],how="left", on="STRUCTURE_NAME")
            df['UnitCell'] = df['UnitCell'].fillna("1 1 1")
            if Threshold_volume > 0:
                df['Volume [nm^3]'] = df['Volume [nm^3]'].fillna(Threshold_volume)
                mask = df['Volume [nm^3]'] <= Threshold_volume
                rejected_threshold = df[~mask]['STRUCTURE_NAME'].tolist()
                if len(rejected_threshold) > 0:
                    print('/!\ The following structures were rejected because they were above volume threshold:', rejected_threshold)
                    df = df[mask]
            if type_ in self.SIMULATION_TYPES['RASPA2']+self.SIMULATION_TYPES['CPP']:
                self.data = df[['STRUCTURE_NAME','UnitCell']].to_records(index=False)
            elif type_ in self.SIMULATION_TYPES['ZEO++']:
                df['ProbeRadius'] = probe_radius
                df['ProbeRadius'] = df['ProbeRadius'].astype(str)
                self.data = df[['STRUCTURE_NAME','ProbeRadius']].to_records(index=False)
            elif type_ in self.SIMULATION_TYPES['HOME']:
                self.atoms = '|'.join([mol2atoms[molecule] for molecule in MOLECULES])
                self.cutoff = float(cutoff)
                df['supercell_wrap'] = df['UnitCell'].apply(lambda x: x.replace(' ','|'))
                self.data = df[['STRUCTURE_NAME','supercell_wrap']].to_records(index=False)
                self.home = True

        pd.DataFrame({'Structures':[],"CPU_time (s)":[]}).to_csv(os.path.join(self.OUTPUT_PATH,"time.csv"), index=False)
        print("%s simulation of %s"%(type_,' '.join(MOLECULES)))


    def generate_files(self, path_to_work, type_, mol2atoms, molecule_dict={}, **kwargs):
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

            path_to_FRAMEWORK = os.path.join(MATSCREEN, "Raspa_screening_templates/FRAMEWORK_%s"%type_)
            if not os.path.isfile(path_to_FRAMEWORK):
                path_to_FRAMEWORK = os.path.join(MATSCREEN, "Raspa_screening_templates/FRAMEWORK")
            if type_ == 'nvt':
                self.pressures = ["0"]
            FRAMEWORK_file = self.generate(path_to_FRAMEWORK, **kwargs)
            if self.movie:
                FRAMEWORK_file += "\nMovies yes\nWriteMoviesEvery "+str(kwargs['N_print'])

            for (index, (key, value)) in enumerate(molecule_dict.items()):
                if len(mol2atoms[key].split(' ')) == 1:
                    rotproba = "0.0"
                elif value is True:
                    rotproba = "0.2"
                else:
                    rotproba = "1.0"
                component = dedent("""
                Component %d MoleculeName                     %s
                            MoleculeDefinition               TraPPE
                            RotationProbability              %s
                """%(index, key, rotproba))
                if type_ == 'coad':
                    component += dedent("""
                                TranslationProbability           0.5
                                IdentityChangeProbability        1.0
                                NumberOfIdentityChanges          2
                                IdentityChangesList              0 1
                                SwapProbability                  1.0
                                CreateNumberOfMolecules          0
                                MolFraction                      %s
                    """%value)
                elif type_ == 'ads' or type_ == 'samplerestart':
                    component += "            CreateNumberOfMolecules"+(" 1" if type_ == 'samplerestart' else " 0")+'\n'
                    if value:
                        component += "            ExtraFrameworkMolecule           yes\n"
                        component += "            TranslationProbability           0.2\n"
                    else:
                        component += "            ExtraFrameworkMolecule           no\n"
                        component += "            TranslationProbability           1.0\n"
                        component += "            SwapProbability                  1.0\n"
                        component += "            BlockPockets                     yes\n"
                        component += "            BlockPocketsFileName             BLOCKPOCKET\n"
                elif type_ == 'pt':
                    component += "            ExtraFrameworkMolecule           yes\n"
                    component += "            TranslationProbability           1.0\n"
                    component += "            ParallelTemperingProbability     0.05\n"
                    component += "            CreateNumberOfMolecules"+(" NUMCATION")*len(self.temperatures)+'\n'
                INPUT_file += component

            if type_ == 'pt':
                for k, temp in enumerate(self.temperatures):
                    INPUT_file += "\nFramework %s\n"%(k) + FRAMEWORK_file + "\nExternalTemperature %s\n"%(temp)
            elif not (type_ == 'sa' or type_ == 'sp' or type_ == 'info'):
                INPUT_file += "\nFramework 0\n" + FRAMEWORK_file + "\nExternalTemperature TEMPERATURE\n"
            self.write_file(INPUT_file, os.path.join(path_to_work, "INPUT"))
            if type_ == 'sp' or type_ == 'sample':
                if not os.path.exists(os.path.join(path_to_work, "Coordinates")):
                    raise FileNotFoundError("Check that Coordinates/ directory exists and it is loaded with results from Voronoi simulations")
                if not os.path.exists(os.path.join(path_to_work, "RestartInitial/System_0")):
                    os.mkdir(os.path.join(path_to_work, "RestartInitial"))
                    os.mkdir(os.path.join(path_to_work, "RestartInitial/System_0"))
                RUN_file = open(os.path.join(MATSCREEN, "Raspa_screening_templates/run_%s"%type_), "r").read()
            elif type_ == "diffusion":
                RUN_file = self.generate(os.path.join(MATSCREEN, "Raspa_screening_templates/run_diffusion"), **kwargs)
                os.system("mkdir %s/Output"%path_to_work)
            elif type_ == 'pt':
                RUN_file = open(os.path.join(MATSCREEN, "Raspa_screening_templates/run_pt"), "r").read()
            else:
                RUN_file = open(os.path.join(MATSCREEN, "Raspa_screening_templates/run"), "r").read()

            self.path_to_run = os.path.join(path_to_work, "run")
            self.write_file(RUN_file, self.path_to_run)
            os.system('chmod +x %s'%(self.path_to_run))
            if type_ != 'grid' and type_ != 'samplerestart':
                julia_DATA_path = os.path.join(MATSCREEN, "Raspa_screening_templates/data_%s.jl"%type_)
                if os.path.exists(julia_DATA_path) and not os.system("julia -v"):
                    DATA_path = julia_DATA_path
                    DATA_ext = "jl"
                else:
                    DATA_path = os.path.join(MATSCREEN, "Raspa_screening_templates/data_%s.sh"%type_)
                    DATA_ext = "sh"
                DATA_file = self.generate(DATA_path, **kwargs)
                self.write_file(DATA_file, os.path.join(path_to_work, "data."+DATA_ext))
                if type_ == 'info':
                    merge_file = self.generate(os.path.join(MATSCREEN, "Raspa_screening_templates/merge_info.py"), **kwargs)
                    self.write_file(merge_file, os.path.join(path_to_work,"merge_info.py"))

            restartinitial = Path(os.path.join(path_to_work, "RestartInitial"))
            if self.RESTART and not (restartinitial / "System_0").exists():
                restartfolder = os.path.join(path_to_work, "Restart")
                if not os.path.exists(restartfolder):
                    raise FileNotFoundError("No RestartInitial nor Restart subfolder at %s, cannot make a restarted job."%path_to_work)
                print("The Restart subfolder at %s is renamed into RestartInitial prior to starting the simulation."%path_to_work)
                os.system('mv %s %s'%(restartfolder, restartinitial))
            if os.path.exists(os.path.join(path_to_work, "RestartInitial/System_0")):
                toexpand = list((restartinitial / "System_0").glob("*_"))
                for fname in toexpand:
                    for (i, temp) in enumerate(self.temperatures):
                        dir = (restartinitial / ("System_%i"%(i if type_ == 'pt' else 0)))
                        dir.mkdir(exist_ok=True)
                        for p in self.pressures:
                            newname = "%s%.6f_%g" % (fname.name, float(temp), float(p))
                            os.system("cp %s %s" % (fname, dir / newname))
                    os.system("rm %s" % fname)

        elif type_ in self.SIMULATION_TYPES["ZEO++"]:
            path_to_Output = os.path.join(path_to_work, 'Output')
            if not os.path.exists(path_to_Output):
                os.mkdir(path_to_Output)
            if type_ == "voronoi":
                if not os.path.exists(os.path.join(path_to_work, 'Coordinates')):
                    os.mkdir(os.path.join(path_to_work, 'Coordinates'))
                os.system("cp %s %s"%(os.path.join(MATSCREEN, "Zeo++_screening_templates/extract_vertex.py"),path_to_work))
            elif type_ in ["channel","sa","strinfo","pore","volpo","oms"]:
                os.system("cp %s %s" % (os.path.join(MATSCREEN, "Zeo++_screening_templates/%s.py"%type_), path_to_work))
                if type_ in ["channel","sa","pore","volpo"]:
                    os.system("cp %s %s" % (os.path.join(MATSCREEN, "Zeo++_screening_templates/rad_uff_298K.rad"), path_to_work))
            RUN_file = self.generate(os.path.join(MATSCREEN, "Zeo++_screening_templates/run_%s"%type_), **kwargs)
            self.path_to_run = os.path.join(path_to_work, "run")
            self.write_file(RUN_file, self.path_to_run)
            os.system('chmod +x %s'%(self.path_to_run))

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
            os.system('chmod +x %s'%(self.path_to_run))
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

    def has_completed_computation(self, struc, unitcell, temperature, pressure_input):
        outputs = Path(self.OUTPUT_PATH) / "Output"
        if not outputs.exists():
            return False
        if temperature == -1:
            system0 = outputs / "System_0"
            if not system0.is_dir():
                return False
            list0 = os.listdir(system0)
            if len(list0) == 0:
                return False
            temperature = list0[0].split('_')[-2]
        pressures = [pressure_input] if pressure_input != -1 else self.pressures
        for pressure in pressures:
            name = "output_" + struc + '_' + unitcell.replace(' ', '.') + "_%.6f_%.6g.data"%(float(temperature), float(pressure))
            notfound = True
            for syst in os.listdir(outputs):
                if len(syst) <= 7 or syst[:7] != "System_" or not syst[7:].isdigit():
                    continue
                fullfile = outputs / syst / name
                if not fullfile.is_file():
                    continue
                if not os.system('tail %s | grep -q "Simulation finished"'%(fullfile)):
                    if os.system('tail %s | grep -q "ENERGY DRIFT"'%(fullfile)):
                        # simulation finished with no energy drift
                        notfound = False
                        break
                # either the simulation did not finish or there was some energy drift
                return False
            if notfound:
                return False # this file could not be found
        return True

    def make_sequential_command(self, struc, unitcell, temperature, pressure_if_parallel):
        restart_string = 'yes' if self.RESTART else 'no'
        subcommand = "bash %s "%(self.path_to_run) + struc + " \"" + unitcell + "\""
        if temperature != -1:
            subcommand += (' ' + str(temperature))
        if self.PressureParallelism == "raspa":
            return subcommand + (" \"" + ' '.join(self.pressures) + "\" %s %i"%(restart_string, self.N_init))
        if self.PressureParallelism == "parallel":
            return subcommand + (" %s %s %i"%(pressure_if_parallel, restart_string, self.N_init))
        command = subcommand + (" %s %s %i"%(self.pressures[0], restart_string, self.N_init))
        assert self.PressureParallelism == "sequence"
        root = Path(self.path_to_run).parent
        restartinitial = root / 'RestartInitial'
        if not restartinitial.exists():
            os.system("mkdir %s && mkdir %s/System_0"%(restartinitial, restartinitial))
            if self.type_ == 'pt':
                for i in range(1, len(self.temperatures)):
                    os.system("mkdir %s/System_%i"%(restartinitial, i))
        restart = root / 'Restart'
        newNinit = self.N_init / 2
        for i in range(1, len(self.pressures)):
            pressure = self.pressures[i]
            if isinstance(unitcell, str):
                newunitcell = unitcell.replace(' ', '.')
            else:
                newunitcell = unitcell.apply(lambda x: x.replace(' ', '.'))
            if self.type_ != "pt":
                if isinstance(temperature, str):
                    newtemperature = ("%.6f"%(float(temperature))).strip()
                else:
                    newtemperature = temperature.apply(lambda t: "%.6f"%(float(t)))
                fileradical = "restart_" + struc + '_' + newunitcell + '_' + newtemperature
                oldfile = fileradical + '_' + ("%5g"%(float(self.pressures[i-1]))).lstrip()
                newfile = fileradical + '_' + ("%5g"%(float(pressure))).lstrip()
                command += " ; cp %s/System_0/"%(restart) + oldfile + " %s/System_0/"%(restartinitial) + newfile + " && " + subcommand + (" %s %s %i"%(pressure, 'yes', newNinit))
            else:
                for (j, T) in enumerate(self.temperatures):
                    newtemperature = ("%.6f"%(float(T))).strip()
                    fileradical = "restart_" + struc + '_' + newunitcell + '_' + newtemperature
                    oldfile = fileradical + '_' + ("%5g"%(float(self.pressures[i-1]))).lstrip()
                    newfile = fileradical + '_' + ("%5g"%(float(pressure))).lstrip()
                    command += " ; cp %s/System_%i/"%(restart, j) + oldfile + " %s/System_%i/"%(restartinitial, j) + newfile
                command += " && " + subcommand + (" %s %s %i"%(pressure, 'yes', newNinit))
        return command

    def run(self, inputs):
        """A module to run one simulation on a given node according to the current process id
        """

        t0 = time()
        FRAMEWORK_NAME,UNITCELL,TEMPERATURE,PRESSURE = inputs
        command = self.make_sequential_command(FRAMEWORK_NAME, UNITCELL, TEMPERATURE if self.type_ != 'pt' else -1, PRESSURE)
        if len(self.NODES) != 0:
            worker = int(mp.current_process()._identity[0])
            nnode = len(self.NODES)
            index = (worker-1)%nnode
            HOST = self.NODES[index]
            command = "ssh %s '%s'"%(HOST,command)
        os.system(command)
        output_dict = {'Structures':[FRAMEWORK_NAME], "CPU_time (s)":[int(time()-t0)]}
        pd.DataFrame(output_dict).to_csv(os.path.join(self.OUTPUT_PATH,"time.csv"),mode="a",index=False,header=False)


    def run_home(self, inputs):
        """Runs calculations using a homemade python algorithm based on pymatgen
        """

        t0 = time()
        structure_name, supercell_wrap, *_ = inputs
        if self.type_ == "surface_sample":
            supercell_wrap = self.n_sample
        if len(self.NODES) == 0:
            command = "%s %s %s %s %s %s \"%s\" \"%s\" "%(sys.executable, self.path_to_run, structure_name, self.cutoff, self.forcefield, self.temperatures[0], self.atoms, supercell_wrap)
            print(command)
        else:
            worker = int(mp.current_process()._identity[0])
            nnode = len(self.NODES)
            index = (worker-1)%nnode
            HOST = self.NODES[index]
            command = "ssh %s \"%s %s \\\"%s\\\" %s %s %s %s \\\"%s\\\" \""%(HOST, sys.executable, self.path_to_run, structure_name, self.cutoff, self.forcefield, self.temperatures[0], self.atoms, supercell_wrap)
        os.system(command)
        output_dict = {'Structures':[structure_name], "CPU_time (s)":[int(time()-t0)]}
        pd.DataFrame(output_dict).to_csv(os.path.join(self.OUTPUT_PATH,"time.csv"),mode="a",index=False,header=False)

    def trim_completed_glost(self, output):
        with open(output) as f:
            lines = f.readlines()
        if len(lines) == 0:
            return 0
        if (not self.SKIPDONE) or lines[0][:4] != "bash":
            return len(lines)
        num = 0
        with open(output, 'w') as f:
            for line in lines:
                splits = shlex.split(line)
                assert splits[0] == "bash"
                struc = splits[2]
                unitcell = splits[3]
                pressure = splits[-3]
                if len(splits) == 8: # bash run struc unitcell T P restart Ni
                    temperature = splits[4]
                else:
                    temperature = -1
                if ' ' in pressure:
                    pressure = -1
                if not self.has_completed_computation(struc, unitcell, temperature, pressure):
                    f.write(line)
                    num += 1
        return num

    def glost_list(self):
        """ Print out glost list for mpirun"""
        df = pd.DataFrame.from_records(self.data)
        struc = df.iloc[:,0]
        opt = df.iloc[:,1]
        output = os.path.join(self.OUTPUT_PATH, "glost.list")
        if os.path.exists(output):
            os.remove(output)
        temperatures = self.temperatures
        if self.type_ == 'pt':
            temperatures = [-1]
        pressures = self.pressures
        if self.PressureParallelism != "parallel":
            pressures = [-1]
        for pressure in pressures:
            for temperature in temperatures:
                if self.home:
                    command = "%s %s "%(sys.executable, self.path_to_run) + struc + " %s %s %s \"%s\" "%(self.cutoff, self.forcefield, temperature, self.atoms) + " \"" + opt + "\""
                else:
                    command = self.make_sequential_command(struc, opt, temperature, pressure)
                command.to_frame().to_csv(output, index=False, header=False, quoting=csv.QUOTE_NONE, quotechar='', mode='a')
        return self.trim_completed_glost(output)
        #return len(pressures)*len(self.temperatures)*(1 if self.type_ == 'pt' else len(struc))

    def slurm_job(self):
        """ Print out a slurm input that can be given to sbatch"""
        num = min(40, (1 + self.glost_list() // 2))
        output = os.path.join(self.OUTPUT_PATH, "input.slurm")
        with open(output, 'w') as f:
            f.write("""#!/bin/bash
#SBATCH --ntasks=%i             # Nombre total de processus MPI
#SBATCH --ntasks-per-node=%i    # Number of MPI tasks per node
#SBATCH --hint=nomultithread    # Disable hyperthreading
#SBATCH --job-name=%s          # Jobname
#SBATCH --output=%%x-%%j.log      # Output file
#SBATCH --error=%%x-%%j.err       # Error file %%x is the jobname, %%j the jobid
#SBATCH --time=20:00:00         # Expected runtime HH:MM:SS (max 100h)
#SBATCH --account=drd@cpu       # To specify cpu accounting: <account> = echo $IDRPROJ

##SBATCH --nodes=1               # Number of Nodes (uncomment to use)
##SBATCH --cpus-per-task=1       # Number of OpenMP threads (uncomment to use)

##SBATCH --partition=<partition>       # To specify partition (see IDRIS web site for more info)
##SBATCH --qos=qos_cpu-dev       # Uncomment for job requiring less than 2 hours
##SBATCH --qos=qos_cpu-t4        # Uncomment for job requiring more than 20h (only one node)

# Manage modules
module purge
module load parallel

# Execute commands
#srun %s/glost_launch glost.list
parallel < glost.list
"""%(num, num, self.type_, os.environ['GLOST_DIR']))


    def mp_run(self):
        """Using multiprocessing, runs in parallel the simulations
        """
        t0 = time()

        run = self.run_home if self.home else self.run
        temperatures = self.temperatures
        if self.type_ == 'pt':
            temperatures = [-1]
        pressures = self.pressures
        if self.PressureParallelism != "parallel":
            pressures = [-1]
        data = [(FRAMEWORK_NAME,UNITCELL,TEMPERATURE,PRESSURE)
                 for (FRAMEWORK_NAME,UNITCELL) in self.data
                 for TEMPERATURE in temperatures
                 for PRESSURE in pressures
                 if not (self.SKIPDONE and self.has_completed_computation(FRAMEWORK_NAME, UNITCELL, TEMPERATURE, PRESSURE))]

        print(data)
        with mp.Pool(processes=self.nprocs) as p:
            p.map(run, data)

        print("SIMULATIONS COMPLETED after %.2f hour(s)"%((time()-t0)/3600))
