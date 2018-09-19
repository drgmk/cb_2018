Getting started
===============

install photodynam
------------------
* Clone photodynam from https://github.com/dfm/photodynam
* Build by typing >make

install rebound
---------------
* >pip3 install rebound --user

modules
-------
* this needs to be done with each login, or saved with >module save and then restored with >module restore
* >module load GCC/6.4.0-2.28 OpenMPI/2.1.1 Python/3.6.3 IPython/6.2.1-Python-3.6.3 matplotlib astropy

run jupyter
-----------
* >jupyter notebook

Files
=====

funcs.py
--------
Various handy functions

known_system_params.py
----------------------
List of orbital parameter setups for different circumbinary planet systems

model_light_curve.ipynb
-----------------------
Notebook to generate a model light curve using photodynam

planet_search.ipynb
-------------------
Example search for Kepler 16b