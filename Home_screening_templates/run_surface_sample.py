import os
import sys
from time import time
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

accessible_mean_energy, min_energy, boltz_energy = lj_sampler.evaluate_from_surface(structure_name)

pd.DataFrame({"Structure_name":N_atoms*[structure_name], "Adsorbent_name":ATOMS, "Acessible_average_energy":accessible_mean_energy, "Minimum_energy":min_energy, "Boltzmann_average_energy":boltz_energy}).to_csv("home_output.csv", index=False, header=False, mode='a')

print(structure_name, ATOMS, accessible_mean_energy, min_energy, boltz_energy)
