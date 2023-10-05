# %%
import pandas as pd
import numpy as np 

import os
os.system("head -1 -q Output/*.vol > results_vol.data")

df = pd.read_csv('results_vol.data', delim_whitespace=True, header=None)

# %%
df["Structures"] = df.iloc[:,1].str.replace("Output_vol/results_","",regex=True).str.replace(".vol","",regex=True)
df["Density"] = df.iloc[:,5] # g/mL
df["cell_volume"] = df.iloc[:,3] # A^3 ?
df["AV_VF"] = df.iloc[:,9]
df["NAV_VF"] = df.iloc[:,15]
df["VF"] = df["AV_VF"]+df["NAV_VF"]

# %%
df[["Structures", "Density", "cell_volume", "AV_VF", "NAV_VF", "VF"]].to_csv('Volume.csv', index=False)
# %%
