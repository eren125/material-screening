import sys
import pandas as pd 
import numpy as np 

struc = sys.argv[1]
DIST_THRESHOLD = 2.0

def extract_vertex(structure):
    f = open("Output/%s.nt2"%structure,'r')
    lines = f.readlines()
    s = 0
    for i in range(len(lines)):
        if lines[i] == "\n":
            s = i
    lines = lines[1:s]
    df = pd.DataFrame([l.split()[1:9] for l in lines])
    df.rename(columns={0:'x',1:'y',2:'z',3:'dist_to_nearest',4:'atom_1',5:'atom_2',6:'atom_3',7:'atom_4'}, inplace=True)
    return df

df = extract_vertex(struc)

df_out = df[df['dist_to_nearest'].astype(float)>DIST_THRESHOLD].sort_values(by=['dist_to_nearest'])

df_out.to_csv("Coordinates/"+struc+'.csv', index=False, sep=" ", encoding='utf-8')
