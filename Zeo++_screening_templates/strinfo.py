# %%
import pandas as pd
import numpy as np
import chemparse
from collections.abc import Mapping

_MAXCOL = 25

import os
os.system('cat Output/*.strinfo > results_strinfo.data')

# %%
# https://rednafi.github.io/reflections/safer-operatoritemgetter-in-python.html
class safe_itemgetter:
    """
    Return a callable object that fetches the given item(s)
    from its operand.
    """

    __slots__ = ("_items", "_call")

    def __init__(self, item, *items, default=0):
        if not items:
            self._items = (item,)

            def func(obj):
                if isinstance(obj, Mapping):
                    return obj.get(item, default)
                if (item > 0 and len(obj) <= item) or (
                    item < 0 and len(obj) < abs(item)
                ):
                    return default
                return obj[item]

            self._call = func
        else:
            self._items = items = (item,) + items

            def func(obj):
                if isinstance(obj, Mapping):
                    get = obj.get  # Reduce attibute search call.
                    return tuple(get(i, default) for i in items)

                return tuple(
                    default
                    if (i > 0 and len(obj) <= i)
                    or (i < 0 and len(obj) < abs(i))
                    else obj[i]
                    for i in items
                )

            self._call = func

    # ----------------- same as operator.itemgetter --------------#

    def __call__(self, obj):
        return self._call(obj)

    def __repr__(self):
        return "%s.%s(%s)" % (
            self.__class__.__module__,
            self.__class__.__name__,
            ", ".join(map(repr, self._items)),
        )

    def __reduce__(self):
        return self.__class__, self._items

organic = ["C", "H", "O", "N"]
halogen = ["Br", "Cl", "F", "I"]
metalloid = ["As", "B", "Ge", "Te", "Sb", "Si"]
ametal = ["Se", "S", "P"]
noble = ["He", "Ne", "Ar", "Kr", "Xe"]
metal_stable = [
"Li", "Be", 
"Na", "Mg", "Al",
"K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co",
"Ni", "Cu", "Zn", "Ga",
"Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Ru", "Rh",
"Pd", "Ag", "Cd", "In", "Sn",
"Cs", "Ba", "La", "Ce", "Pr", "Nd", "Sm", "Eu",
"Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
"Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
"Tl", "Pb", "Bi", 
]
# https://www.thoughtco.com/list-of-radioactive-elements-608644
radioactive = [
"Tc", 
"Pm",
"Po", "At", "Rn",
"Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am",
"Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
"Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
"Nh", "Fl", "Mc", "Lv", "Ts", "Og"
]

# %%
df = pd.read_csv('results_strinfo.data', delim_whitespace=True, header=None, names=np.arange(0,_MAXCOL))
# %%
Ninit_col = 18
max_frameworks = df[4].max()
if max_frameworks + Ninit_col > _MAXCOL:
    print('Not enough max_col')

# %%
df["Structures"] = df.iloc[:,0].str.replace("Output_info/results_","",regex=True).str.replace(".strinfo","",regex=True)
df["chem_comp"] = df.iloc[:,1]
df['solvent_removed'] = (df[2]==df[4]).astype(int)
df["framework_count"] = df.iloc[:,4].astype(int)
df["framework_mean_dim"] = (df[np.arange(Ninit_col,max_frameworks+Ninit_col)].fillna(0).sum(axis=1)/df["framework_count"]).fillna(-1)

# %%
df[['Structures','chem_comp', 'solvent_removed', 'framework_count', 'framework_mean_dim']].to_csv('zeoinfo.csv', index=False)

# %%
chemcomp_dict = df['chem_comp'].apply(lambda x: chemparse.parse_formula(x))
df['atoms_count'] = chemcomp_dict.apply(lambda x: sum(x.values()))
inv_tot_percent = 100/(df['atoms_count']).to_numpy()

# %%
for symbol in organic:
    df["%s%%"%symbol] = chemcomp_dict.apply(lambda dict: dict.get(symbol,0)).to_numpy() * inv_tot_percent

# %%
df["halogen%"] = chemcomp_dict.apply(lambda dict: sum(safe_itemgetter(*tuple(halogen))(dict))).to_numpy() * inv_tot_percent

df["metalloid%"] = chemcomp_dict.apply(lambda dict: sum(safe_itemgetter(*tuple(metalloid))(dict))).to_numpy() * inv_tot_percent

df["ametal%"] = chemcomp_dict.apply(lambda dict: sum(safe_itemgetter(*tuple(ametal))(dict))).to_numpy() * inv_tot_percent

df["metal%"] = chemcomp_dict.apply(lambda dict: sum(safe_itemgetter(*tuple(metal_stable))(dict))).to_numpy() * inv_tot_percent

# %%
df["M/C"] = df["metal%"]/df["C%"]
df["M/O"] = df["metal%"]/df["O%"]

# %%
df["radioactive%"] = chemcomp_dict.apply(lambda dict: sum(safe_itemgetter(*tuple(radioactive))(dict))).to_numpy() * inv_tot_percent

# %%
descriptors = ['Structures', 'solvent_removed', 'framework_count', 'framework_mean_dim', 'atoms_count', 'C%', 'H%', 'O%', 'N%', 'halogen%', 'metalloid%', 'ametal%', 'metal%', "radioactive%", 'M/C', 'M/O',]

df[descriptors].to_csv('zeoinfo_descriptors.csv', index=False)

# %%
