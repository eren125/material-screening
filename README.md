# Material Screening

A python-based screening procedure for both Raspa2 and Zeo++ softwares

## Installations

This code relies on 3 softwares that need a prior installation step

### Python3.8.5

Download and install Python3.8.5 at https://www.python.org/
Other versions of python3 should be ok, but I have not tested them yet.

### Raspa2

Follow the instructions of the **Raspa2 Github repository** to compile it:
https://github.com/iRASPA/RASPA2
Remember the path to the directory containing the ```bin/``` directory containing the binary file ```simulate```

### Zeo++

Download and compile Zeo++ at http://www.zeoplusplus.org/download.html
Remember the path to the directory containing the compiled ```network``` binary file

## Initial set-up

This part explains how you will need to set-up your working environment

### Bash environment

In the ```set_environment``` file, replace the example paths by the accurate ones.

## RASPA environment

All molecules will be taken from ${RASPA_DIR}/share/raspa/molecules/TraPPE, which must
contain at least the definition for helium (as well as all other molecules you need).

If the folder is missing, simply create it.
The file for helium can be taken from https://github.com/numat/RASPA2/blob/master/molecules/TraPPE/helium.def.

Any used force field must be present in ${RASPA_DIR}/share/raspa/forcefields. In particular,
the default force field UFF should be copied from ${MATSCREEN}/forcefield if not already
present.

### Python environment

I used a virtual environment for ```python3.8.5``` but it is optional.
To set-up a virtual environment and install the required packages in ```requirements.txt```, follow those instructions:
```
./path_to_python3.8.5/bin/python3.* -m venv path_to_venv
source path_to_venv/bin/activate
pip install -r requirements.txt
```

If you want to go back to your global python environment type: ``deactivate``

## Database of materials and force fields

### Force Field

For the forcefield, you can use the default ones in Raspa2 add your own force fields adapted to your applications: You can find in the ```forcefield``` directory an exmaple of UFF forcefield built on this <a href="https://pubs.acs.org/doi/10.1021/ja00051a040">paper</a>. You can also build a Dreiding forcefield from this paper <a href="https://pubs.acs.org/doi/10.1021/j100389a010">paper</a>, usually it is mixed with UFF to complete the missing atoms.

### Material Databases

<ul>
    <li><b>CoRE MOF Database</b>: The Computation-Ready, Experimental Metal--Organic Framework (MOF) database is an experimental based dabase described in this <a href="https://pubs.acs.org/doi/abs/10.1021/acs.jced.9b00835">paper</a> and can be downloaded via this <a href="https://zenodo.org/record/3677685#.YFnEu3VKhhE">link</a> (about 14,000 structures)</li>
    <li><b>QMOF Database</b>: Quantum MOF (QMOF) database is a DFT-optimized MOFs database <a href="https://chemrxiv.org/articles/preprint/Machine_Learning_the_Quantum-Chemical_Properties_of_Metal_Organic_Frameworks_for_Accelerated_Materials_Discovery_with_a_New_Electronic_Structure_Database/13147616/1?file=25304507">paper</a> and can be downloaded via this <a href="https://figshare.com/articles/dataset/QMOF_Database/13147324">link</a> (about 14,000 structures)</li>
    <li><b>MOFDB</b>: ToBaCCo and hMOF datasets can be downloaded via this <a href="https://mof.tech.northwestern.edu/databases">link</a>. The Topologically-Based Crystal Constructor generates automatically MOFs based on the metal, the linker and the topology. The algorithm is described in detail in this <a href="https://pubs.acs.org/doi/abs/10.1021/acs.cgd.7b00848">paper</a> and the hMOF database is another hypothetical MOF database described in this <a href="https://www.nature.com/articles/nchem.1192">paper</a> </li>
</ul>

Zeolite, Porous Polymer Networks and Covalent Organic Framework can also be found as databases.

## Running your first simulations

Using `screen.py` command line tool to generate adsorption or coadsorption simulation files and then run them.

Once the structures have been decided and put in `${RASPA_DIR}/share/raspa/structures/cif/`,
put their names in a one-column CSV file with "Structures" as its only header.
Then, execute `screen.py` with the `-t info` option to extract the information necessary
to run subsequent simulations. You should then launch `./data.sh` to obtain an `info.csv`
file which should replace the one in `${MATSCREEN}/data/`.

Example of usage in `job_example.sh`

Need help? Execute: `screen.py --help`.
