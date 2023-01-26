# %%
import pandas as pd

import os
os.system("cat Output/*.oms > results_oms.data")

df = pd.read_csv('results_oms.data', delim_whitespace=True, header=None)

# %%
df['Structures'] = df.iloc[:,0].str.replace("Output/results_","",regex=True).str.replace(".oms","",regex=True)
df.rename(columns={2:"oms_count"}, inplace=True)
# OMS density could be a good descriptor 

# %%
df[['Structures','oms_count']].to_csv('oms.csv', index=False)

# %%
