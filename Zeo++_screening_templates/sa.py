# %%
import pandas as pd
import numpy as np 

import os
os.system("head -1 -q Output/*.sa > results_sa.data")

df = pd.read_csv('results_sa.data', delim_whitespace=True, header=None)

# %%
df["Structures"] = df.iloc[:,1].str.replace("Output_sa/results_","",regex=True).str.replace(".sa","",regex=True)
df["ASA_m2/cm3"] = df.iloc[:,9]
df["NASA_m2/cm3"] = df.iloc[:,15]
df["SA_m2/cm3"] = df["ASA_m2/cm3"]+df["NASA_m2/cm3"]

# %%
df[["Structures", "ASA_m2/cm3", "NASA_m2/cm3", "SA_m2/cm3"]].to_csv('SurfaceArea.csv', index=False)
# %%
