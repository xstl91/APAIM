*APAIM* - lateral interaction calculation
========================


Basic information
------------------------

Here we provide a relatively universal code to implement the **Augmented Pairwise Additive Interaction Model** we've proposed in the article:

> Augmented pairwise additive interaction model for lateral adsorbate interactions: the NO-CO reaction system on Rh(100) and Rh(111)

The code is developed based on Python 2.7.  
No external module is needed unless one wants to visualize the configurations, in which case *numpy* and *matplotlib* modules are needed.

The code is restricted to:
* (100) and (111) facets of FCC crystal, or other surfaces with similar structure
* 4×4 unit cell, though the model itself can be extended


Contains
-------------------------

* `APAIM`: a module to implement the **Augmented Pairwise Additive Interaction Model**
* `100_site.tif` and `111_site.tif`: the notation used for surface sites
* `dataset`: testing datasets for the NO-CO reaction system on Rh(100) and Rh(111)
* `README.md`: a brief introduction
* Some usage examples:
  * `single.py`: get the predicted binding energy for a configuration
  * `dataset_error.py`: analyze the error information for a dataset
  * `dataset_cutoff.py`: analyze the interaction-cutoff-distance information for a dataset
  * `stable.py`: get a local minimum using the model Hamiltonian built


*APAIM* module
-------------------------

#### Usage
Use the `help` function in Python to know more about the usage of this module:
```python
$ python
/* version information */
>>> from APAIM import *       # import the module
/* choose a surface */
>>> help('APAIM')             # get usage information of the module
>>> help(Config)             # get usage information of the class developed for single configuration
>>> help(function_name)      # get usage information of the functions developed
```

#### Current system
The APAIM module is currently set to the NO-CO reaction system on Rh(100) and Rh(111).  
Isolated binding energies are provided in directory `APAIM/p_single`.  
Interaction energies are provided in directory `APAIM/p_inter`.  
All energies provided here are in eV.

#### Extension
To reuse the module to other surface adsorption systems:
* modify the following constants in `APAIM/constant.py`:
  * `species_all` for adsorbate types
  * `recon_pair_hh1` and `recon_pair_hh2` for reconstruction on (100) surface
* provide corresponding isolated binding energies
* provide corresponding interaction energies 
