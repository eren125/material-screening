# Raspa_filegen

raspa file generator using python

## Installation of Raspa

Follow the instructions of the **Raspa2 Github repository**:
https://github.com/iRASPA/RASPA2

## Initial set-up

Set the Raspa environment variables using the following command in this directory:
```
source set_environment
```
In order not to run it every time, you can put the command ```source path_to_dir\set_environment``` in you bashrc file ```vi ~/.bashrc``` on linux

### Python environment

I used a virtual environment of ```python3.8.5``` with the packages of requirements.txt
```
./path_to_python3.8.5/bin/python3.* -m venv path_to_venv
source path_to_venv/bin/activate
pip install -r requirements.txt
```

If you want to go back to your global environment type: deactivate
If you want the virtual environment to activate automatically when you open the terminal copy and paste ```source path_to_venv/bin/activate``` into ```~/.bashrc``` (for linux) ```~/.profile``` (for MacOs) and for Windows I suggest you to install WSL (Linux Sub-system for Windows).

### Raspa2 environment

Change ```set_environment``` file so that it corresponds to the Path of your own Raspa ```Simulations/``` directory.
Replace ```/home/emmanuel/Simulations``` by your own path.

## Running your first simulations

Using ```raspagen.py``` command line tool to generate adsorption or coadsorption simulation files and then run them

Example of usage in ```run-python.sh```
