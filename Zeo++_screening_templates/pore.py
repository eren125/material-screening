# %%
import pandas as pd

import os
os.system("cat Output/*.res > results_pore.data")

df = pd.read_csv('results_pore.data', delim_whitespace=True, header=None)
# %%
df['Structures'] = df[0].str.replace('Output/results_','',regex=True).str.replace(".res", "",regex=True)
df.rename(columns={1:"D_i_vdw_uff298", 2:"D_f_vdw_uff298", 3:"D_if_vdw_uff298"}, inplace=True)

# %%
df[["Structures", "D_i_vdw_uff298", "D_f_vdw_uff298", "D_if_vdw_uff298"]].to_csv('PoreSize.csv', index=False)
# %%
