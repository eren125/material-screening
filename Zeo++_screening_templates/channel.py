# %%
import pandas as pd
import numpy as np 

import os
os.system("head -1 -q Output/*.chan > results_chan.data")

max_col = 40
df = pd.read_csv('results_chan.data', delim_whitespace=True, header=None , names=np.arange(0,max_col))

# %%
Ncol_init = 6
max_chan = df[1].max()
if max_chan+Ncol_init > max_col:
    print('Change max column stting to a higer number')

# %%
df["Structures"] = df.iloc[:,0].str.replace("Output/results_","",regex=True).str.replace(".chan","",regex=True)
df["chan_count"] = df.iloc[:,1]
df["chan_mean_dim"] = (df[np.arange(Ncol_init,max_chan+Ncol_init)].fillna(0).sum(axis=1)/df["chan_count"]).fillna(-1)

# %%
df[['Structures', 'chan_count', 'chan_mean_dim']].to_csv('Channel.csv', index=False)
# %%
