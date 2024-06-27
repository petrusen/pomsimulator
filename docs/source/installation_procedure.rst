Installation procedure
======================

POMSimulator can be installed via git and pip (see below).

.. code:: bash

     $ git clone <pomsimulator-repo>
     $ cd pomsimulator
     $ git checkout v1.0.0 # stable version - release 1.0
     $ pip install -e .

.. attention:: Installing POMSimulator directly on your /home/user may result in incompatibilities with other Python dependent packages. Thus, we highly recommend creating a separate Python environment.

.. code:: bash

     $ python3 -m venv pomsim_env
     $ cd pomsim_env
     $ source ./bin/activate
     (pomsim_env)$ git clone <pomsimulator-repo>
     (pomsim_env)$ cd pomsimulator
     (pomsim_env)$ git checkout v1.0.0 # stable version - release 1.0
     (pomsim_env)$ pip install -e .   #it will automatically install the required dependencies

To check whether the installation was successful, we recommend running a short simulation (i.e., consisting of monomers, dimers and trimers) already prepared in the package:

.. code:: bash

     (pomsim_env)$  python3 pomsimulator/simulations/simulation_test.py


.. note:: Installation of POMSimulator via PyPi is underway
