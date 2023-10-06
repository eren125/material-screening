import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

D_min = 3
i_min = D_min*10
def import_psd(path_file, plot=True):
    Lines = open(path_file, "r").readlines()
    if not Lines:
        return [], []
    Content = np.array([line.split() for line in Lines[12:]])
    X = Content[i_min:200,0].astype(float)
    Y = Content[i_min:200,1].astype(int)
    Xs = X[9::5]
    if plot:
        fig = plt.figure()
        ax = fig.add_axes([0,0,1,1])
        ax.bar(X,Y,0.1)
        plt.xticks(Xs, Xs, rotation='vertical')
        plt.xlim(left=3)
    return X,Y


def pore_dist_moments(X,Y,threshold=0.1):
    if len(X)==0 or len(Y)==0:
        return np.nan, np.nan, np.nan, np.nan
    N = sum(Y)
    mean = (X*Y).sum()/N
    std = np.sqrt(((X-mean)**2*Y).sum()/N)
    skewness = (((X-mean)/std)**3*Y).sum()/N
    kurtosis = (((X-mean)/std)**4*Y).sum()/N
    return mean, std, skewness, kurtosis

L_1, L_2, L_3, L_4 = [], [], [], []
L_struc = []
directory = "Output"
for filename in os.listdir(directory):
    X, Y = import_psd(os.path.join(directory,filename), plot=False)
    mean, std, skewness, kurtosis = pore_dist_moments(X,Y,threshold=0.05)
    L_1.append(mean)
    L_2.append(std)
    L_3.append(skewness)
    L_4.append(kurtosis)
    L_struc.append(filename.replace("results_","").replace(".psd",""))

df_psd = pd.DataFrame(data={"Structures":L_struc, "pore_dist_mean":L_1, "pore_dist_std":L_2, "pore_dist_skewness":L_3, "pore_dist_kurtosis":L_4})

columns = ["Structures","pore_dist_mean", "pore_dist_std", "pore_dist_skewness", "pore_dist_kurtosis"]

df_psd[columns].to_csv("psd_%s.csv"%D_min,index=False)
