[![DOI](https://img.shields.io/badge/DOI-10.1039/D0SC03530K-orange)](https://doi.org/10.1039/D0SC03530K)
[![Latest Version](https://img.shields.io/github/v/release/petrusen/pomsimulator)](https://github.com/petrusen/pomsimulator/releases/latest)
[![License](https://img.shields.io/badge/license-AGPL3.0-darkgreen)](https://github.com/petrusen/pomsimulator/blob/main/LICENSE.txt)

![](docs/.img/pomsimulator_logo.png)

# POMSimulator

- [Introduction](#Introduction)
- [Installation](#Installation)
- [Quickstart](#Quickstart)
- [How to Cite](#Howtocite)
- [License and Copyright Information](#licenseandcopyrightinformation) 
- [Support and Contact](#supportandcontact)

## Introduction 

POMSimulator is a software for predicting the aqueous speciation and self-assembly mechanism of polyoxometalates. Based on a pure Python framework, the method 
generates reaction maps using Graph Theory, and solves the non-linear equations related to the speciation models. The software creates a collection of formation constants for
all the species in the simulation, and a list of the chemical reactions. The package also contains a set of python scripts for analysing and
plotting this data.

**Authors**: Enric Petrus, Mireia Segado-Centellas, Carles Bo

**Institution**: research group of Carles Bo, at the Institute of Chemical Research of Catalonia (ICIQ) - Tarragona, Spain

**Developer team**: Jordi Buils, Diego Garay-Ruiz, Enric Petrus

## Installation

POMSimulator can be installed via git and pip (see below).

```console
git clone <pomsimulator-repo>
cd pomsimulator
pip install -e .
```

To ensure that the installation has been successful, we recommend to run the test `python simulation_test.py` from inside the Simulations folder. 

## Quickstart 

The main simulation file of POMSimulator is stored in `./Simulations`. The code reads the quantum mechanical outputs in `./Inputs` and generates two
output files in the folder `./Outputs`. These two files correspond to the chemical reaction network, and the list predicted of formation constants. Hitherto, the code parses the input data from the calculations of the SCM ADF 2019 package. Nonetheless, the input data can also be introduced separately. 

The run the code, the user has to set the simulation.py file, and the execute it. For example: `python simulation_tungstates.py`. 

With the obtained results, the package offers a set of tools, in the folder `./Utilities`, for analysing and plotting the results.
Before plotting the speciation and/or phase diagrams, it is necessary to linearly scale the predicted formation constants. This is done with 
`scale_constants.py`. The scaling parameters are stored as a text file (`scaling_params.pomsim`), and then employed by `speciation_diagram.py`.

![](docs/.img/pom_workflow.png)

The user can find more information about the basis of the method in the PhD thesis of the author. The document is available in the [TDX-repository](https://www.tesisenred.net/handle/10803/687274). The user can 
also find a written manual in `.docs/manual.md`


## How to Cite

When publishing results generated with POMSimulator, we kindly request you to cite the original paper, and if possible the corresponding release as archived on [Zenodo](https://zenodo.org/records/10689769):

> Petrus, E.; Segado, M.; Bo, C. *Chem. Sci.*, **2020**, 11, 8448-8456 [doi.org/10.1039/D0SC03530K](https://doi.org/10.1039/D0SC03530K)

List of peer-reviewed articles featuring POMSimulator:
* Petrus, E.; Bo, C. *J. Phys. Chem. A.* **2021**, 125, 23, 5212-5219 [doi.org/10.1021/acs.jpca.1c03292](https://doi.org/10.1021/acs.jpca.1c03292)
* Petrus, E.; Segado, M.; Bo, C. *Inorg. Chem.* **2022**, 61, 35, 13708-13718 [doi.org/10.1021/acs.inorgchem.2c00925](https://pubs.acs.org/doi/abs/10.1021/acs.inorgchem.2c00925)
* Petrus, E.; Garay-Ruiz, D.; Reiher, M.; Bo, C. *J. Am. Chem. Soc*, **2023**, 145, 34, 18920-18930 [doi.org/10.1021/jacs.3c05514](https://pubs.acs.org/doi/full/10.1021/jacs.3c05514)
* Garay-Ruiz, D.; Petrus, E.; Bo, C. *Revista de la Societat Catalana de Quimica* **2023**, 22, 23-38 [doi.org/10.2436/20.2003.01.142](https://revistes.iec.cat/index.php/RSCQ/article/view/150830/148565) 
* Petrus, E.; Buils, J.; Garay-Ruiz, D.; Segado-Centellas, M.; Bo, C. *J. Comput. Chem.* **2024**, 1-9 [doi.org/10.1002/jcc.27389](https://doi.org/10.1002/jcc.27389) 

## License and Copyright Information

POMSimulator is distributed under the GNU AFFERO GENERAL PUBLIC LICENSE. For more license and copyright information, see the file `LICENSE.txt`

## Support and Contact

In case you should encounter problems or bugs, please write a short message to one of the following addresses:
[enricpz@icloud.com](enricpz@icloud.com), [jbuils@iciq.es](jbuils@iciq.es), [dgaray@iciq.es](dgaray@iciq.es).
