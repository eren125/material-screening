import pandas as pd

cutoff = 12

lines = open('boxlengths.txt').readlines()
df_box = pd.DataFrame([line.split() for line in lines])
df_box['STRUCTURE_NAME'] = df_box[0].str.replace('Output/System_0/output_','').str.replace('_1.1.1_298.000000_0.data-','')

df_box.rename(columns={1:'x_box [A]',2:'y_box [A]',3:'z_box [A]'},inplace=True)
df_box['a'] = ((cutoff*2)/df_box['x_box [A]'].astype(float)).astype(int) + 1
df_box['b'] = ((cutoff*2)/df_box['y_box [A]'].astype(float)).astype(int) + 1
df_box['c'] = ((cutoff*2)/df_box['z_box [A]'].astype(float)).astype(int) + 1

df_box['UnitCell'] = df_box['a'].astype(str) + ' ' + df_box['b'].astype(str) + ' ' + df_box['c'].astype(str)

df_box.drop(columns=[0,'a','b','c'],inplace=True)

lines = open('volume.txt').readlines()
df_vol = pd.DataFrame([line.split() for line in lines])
df_vol.rename(columns={0:"STRUCTURE_NAME"},inplace=True)
df_vol['STRUCTURE_NAME'] = df_vol['STRUCTURE_NAME'].str.replace('Output/System_0/output_','').str.replace('_1.1.1_298.000000_0.data:volume','')
df_vol.drop(columns=[1,2,3,5],inplace=True)
df_vol.rename(columns={4:'Volume [nm^3]'},inplace=True)
df_vol['Volume [nm^3]'] = df_vol['Volume [nm^3]'].astype(float) / 1000

df_vol_box = pd.merge(df_vol,df_box,how='inner',on='STRUCTURE_NAME')

lines = open('matrix.txt').readlines()
L = [line.split() for line in lines]
S = [x[0] for x in L[0::3]]
M1 = ["%s %s %s"%(x[1],x[2],x[3]) for x in L[0::3]]
M2 = ["%s %s %s"%(x[1],x[2],x[3]) for x in L[1::3]]
M3 = ["%s %s %s"%(x[1],x[2],x[3]) for x in L[2::3]]
df_mat = pd.DataFrame(data={"STRUCTURE_NAME":S, "unit vector a":M1, "unit vector b":M2, "unit vector c":M3})
df_mat['STRUCTURE_NAME'] = df_mat['STRUCTURE_NAME'].str.replace('Output/System_0/output_','').str.replace('_1.1.1_298.000000_0.data-','')
# TODO box lengths can be calculated from the matrix

df = pd.merge(df_vol_box,df_mat,how='inner',on='STRUCTURE_NAME')

lines = open('framework.txt').readlines()
L = [line.split() for line in lines]
S = [x[0] for x in L[0::2]]
M = [x[2] for x in L[0::2]]
D = [x[2] for x in L[1::2]]
df_frame = pd.DataFrame(data={"STRUCTURE_NAME":S, "Framework Mass [g/mol]":M, "Framework Density [kg/m^3]":D})
df_frame['STRUCTURE_NAME'] = df_frame['STRUCTURE_NAME'].str.replace('Output/System_0/output_','').str.replace('_1.1.1_298.000000_0.data:Framework','')

df_merge = pd.merge(df,df_frame,how='inner',on='STRUCTURE_NAME')

df_merge.to_csv('info.csv', index=False, encoding='utf-8')
