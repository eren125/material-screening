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

### Glost

Download and compile Glost at https://github.com/cea-hpc/glost
It is used to run the glost list generated when using the ```-X glost``` of ```screen.py```.

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
    <li><b>MOFXDB</b>: ToBaCCo and hMOF datasets can be downloaded via this <a href="https://mof.tech.northwestern.edu/databases">link</a>. The Topologically-Based Crystal Constructor generates automatically MOFs based on the metal, the linker and the topology. The algorithm is described in detail in this <a href="https://pubs.acs.org/doi/abs/10.1021/acs.cgd.7b00848">paper</a> and the hMOF database is another hypothetical MOF database described in this <a href="https://www.nature.com/articles/nchem.1192">paper</a> </li>
</ul>

Zeolite, Porous Polymer Networks and Covalent Organic Framework can also be found as databases.

## Running your first simulations

Using `screen.py` command line tool to generate adsorption or coadsorption simulation files and then run them.

Once the structures have been decided and put in `${RASPA_DIR}/share/raspa/structures/cif/`,
put their names in a one-column CSV file with "Structures" as its only header.
Then, execute `screen.py` with the `-t info` option to extract the information necessary
to run subsequent simulations. You should then launch `./data.sh` to obtain an `info.csv`
file which should replace the one in `${MATSCREEN}/data/`.

## Workflow example: simulating adsorption in cationic zeolite

Choose your target cationic zeolites and put their CIF files in
`${RASPA_DIR}/share/raspa/structures/cif/`
Store their name in a one-column CSV file with "Structures" as its only header, like in the
previous paragraph. Give that file a name, for instance "foobar.csv".

Then, prepare electrostatic and VdW grids for each structure, to reduce the computation
required subsequently: to do so, run
```
$MATSCREEN/screen.py -ppn X -n X -N 1 -s $MATSCREEN/data/foobar.csv -f FF -x "EXTRA" -m MOLECULES -t grid
```
where:
- `X` is the number of processes to run in parallel.
- `FF` is the name of the force field, which should be placed in `${RASPA_DIR}/share/raspa/forcefield`.
- `EXTRA` is additional command to add in each RASPA INPUT file. For instance,
  `RemoveAtomNumberCodeFromLabel yes` is useful if your CIF file numbers each atom.
  If you do not have any command to add, you can remove the `-x "..."` flag.
- `MOLECULES` is the list of molecules that will be added in the framework. For instance,
  putting `Na CO2 N2` is necessary to simulate coadsorption of CO2 and N2 on zeolites with
  sodium cation. If Na cations are already part of the framework, they should be removed
  from the list.

If the cations are not part of the framework, they should be equilibrated. This can be done
with parallel tempering using the following command:
```
$MATSCREEN/screen.py -ppn X -n X -Ni INIT -N RUN -s $MATSCREEN/data/foobar.csv -f FF -x "EXTRA" -m CATION -t pt -T TEMPERATURES -M
```
where, in addition to the previous replacements:
- `INIT` is the number of initialization steps.
- `RUN` is the number of running steps (for which the movie will be recorded thanks to the `-M` option).
- `CATION` is the name of the cation.
- `TEMPERATURES` is the list of temperatures used for parallel tempering.

Different cation placements can be extracted from the movie, and should be transformed into
restart files for later use. If only the last cation placement is required, you can remove the
`-M` flag and directly take the restart file in the `Restart` folder.

Running `./data.sh` (or `./data.jl` if julia is available) outputs digests of RASPA's
output giving the energy at each recorded step in the Movie, which is useful to check
convergence. These are stored in `DATAenergy...` files.

Once cations have a starting position, the adsorption simulation itself can be run.
- If the cations are fixed, then the most efficient strategy consists in making a new CIF
  file containing the cations
- Otherwise, if cations are mobile, their starting position should be given as a restart
  file to be put in `RestartInitial/System_0`. Either one restart file should be made for
  each appropriate temperature and pressure, or a single restart file can be put whose name
  should be truncated to a trailing underscore, leaving the temperature and pressure blank.
  This file will then be copied and the name expanded for each combination of temperature
  and pressure.

The command to run to perform adsorption simulation is
```
$MATSCREEN/screen.py -ppn X -n X -Ni INIT -N RUN -s $MATSCREEN/data/foobar.csv -f FF -x "EXTRA" -m MOLECULES -c CODE -t ads -T TEMPERATURES -p PRESSURES
```
where
- `CODE` is a list of `0` and `1`, one for each molecule in `MOLECULES`: 1 indicate that
  this molecule is the cation, and 0 indicates that it is a gas whose adsorption is measured.
  If all molecules are gas (the cation is part of the framework), this `-c ...` flag can
  be removed.
- `PRESSURES` is the list of pressures for which the simulation is carried.

Running `./data.sh` (or `./data.jl` if julia is available) outputs a digest of RASPA's
output giving the loading at each recorded step in the Movie. These are stored in
`DATAloading...` files.

In all cases, if a RASPA computation is interrupted, it can be restarted by moving the
`Restart` subfolder into a new folder, renaming the subfolder into `RestartInitial` and
then re-running the same command with an additional `-R` flag in this new folder.
If there is no `RestartInitial` in the current folder and a `-R` command is issued, the
`Restart` folder will be renamed into `RestartInitial`. Beware that previous simulation
contained in the new folder will be overwritten by the new simulation so only restart
a simulation in an occupied folder if you do not need these data.

Additionally, any command can generate files instead of directly running by specifying a
`-X` option:
- `-X exe` is the default option which runs the computation immediately.
- `-X glost` generates a glost input file.
- `-X slurm` generates a slurm input file.

## Other

Example of usage in `job_example.sh`

Need help? Execute: `screen.py --help`.
