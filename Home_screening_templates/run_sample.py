import os
import sys
from time import time, sleep
import numpy as np
import pandas as pd

os.chdir(os.path.dirname(os.path.abspath(__file__)))
SOURCE_DIR = os.environ['MATSCREEN']
sys.path.append(SOURCE_DIR)
from src.ljsampler import load

Args = sys.argv 

ATOMS = Args[1].split('|')
N_atoms = len(ATOMS)
lj_sampler = load(ATOMS, forcefield=Args[2], temperature=float(Args[3]), cutoff=float(Args[4]))

structure_name = Args[5]
supercell = [int(c) for c in Args[6].split('|')]
# lattice_wrap = np.array(Args[7].split('|')).astype(float)
# lattice_matrix = np.array( [supercell[0]*lattice_wrap[:3],supercell[1]*lattice_wrap[3:6],supercell[2]*lattice_wrap[6:9]] )
# example: 5.5673|0|0|0|11.8214|-2.7277|0|0|22.38955

accessible_mean_energy, min_energy, boltz_energy = lj_sampler.evaluate_from_coordinates(structure_name, supercell)

path = os.path.join('.output_written.tmp')
while not os.path.exists(path):
    sleep(0.00001)
os.remove(path)
pd.DataFrame({"Structure_name":N_atoms*[structure_name], "Adsorbent_name":ATOMS, "Acessible_average_energy":accessible_mean_energy, "Minimum_energy":min_energy, "Boltzmann_average_energy":boltz_energy}).to_csv("home_output.csv", index=False, header=False, mode='a')
open(path, 'a')

print(structure_name, ATOMS, accessible_mean_energy, min_energy, boltz_energy)
