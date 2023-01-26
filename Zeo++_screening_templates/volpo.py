# %%
import pandas as pd
import numpy as np 

import os
os.system("head -1 -q Output/*.volpo > results_volpo.data")

df = pd.read_csv('results_volpo.data', delim_whitespace=True, header=None)

# %%
df["Structures"] = df.iloc[:,1].str.replace("Output_volpo/results_","",regex=True).str.replace(".volpo","",regex=True)
df["Density"] = df.iloc[:,5] # g/mL
df["cell_volume"] = df.iloc[:,3] # A^3 ?
df["POA_VF"] = df.iloc[:,9]
df["PONA_VF"] = df.iloc[:,15]
df["PO_VF"] = df["POA_VF"]+df["PONA_VF"]

# %%
df[["Structures", "Density", "cell_volume", "POA_VF", "PONA_VF", "PO_VF"]].to_csv('VolumeOccupiable.csv', index=False)
# %%
