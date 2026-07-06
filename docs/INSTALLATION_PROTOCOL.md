POMSimulator Installation Protocol
* Jordi Buils

July 2026


# Table of Contents

- [Introduction](#introduction)
- [Route 1: Standard Python Installation with venv and pip](#route-1-standard-python-installation-with-venv-and-pip)
  - [Install Python 3.10, 3.11, or 3.12](#install-python-310-311-or-312)
    - [Linux](#linux)
    - [macOS](#macos)
    - [Windows](#windows)
  - [Install Git](#install-git)
    - [Linux](#linux-1)
    - [macOS](#macos-1)
    - [Windows](#windows-1)
  - [Clone POMSimulator and Switch to the gui_alpha Branch](#clone-pomsimulator-and-switch-to-the-gui_alpha-branch)
  - [Create and Activate a Virtual Environment](#create-and-activate-a-virtual-environment)
    - [Linux or macOS](#linux-or-macos)
    - [Windows PowerShell](#windows-powershell)
    - [Windows Command Prompt](#windows-command-prompt)
  - [Install POMSimulator](#install-pomsimulator)
  - [Verify the Installation](#verify-the-installation)
  - [Run the GUI](#run-the-gui)
  - [Minimal Command Summary](#minimal-command-summary)
- [Route 2: Conda or Miniconda Installation](#route-2-conda-or-miniconda-installation)
  - [Install Miniconda](#install-miniconda)
  - [Clone Repository](#clone-repository)
  - [Create and Activate Conda Environment](#create-and-activate-conda-environment)
  - [Install POMSimulator](#install-pomsimulator-1)
  - [Verify Installation](#verify-installation)
  - [Run the GUI](#run-the-gui-1)
  - [Minimal Conda Summary](#minimal-conda-summary)
- [Final Installation Checklist](#final-installation-checklist)
- [Troubleshooting & Error Messages](#troubleshooting--error-messages)


# Introduction

This document explains how to install and run the GUI-enabled version of POMSimulator from the `gui_alpha` branch.

The protocol is divided into two independent installation routes:

1. Standard Python installation using Python 3.10, 3.11, or 3.12,
`venv`, and `pip`.
2. Conda or Miniconda installation using a Conda environment and `pip`
inside that environment.

For both routes, the final GUI launch command is:

python GUI.py


Run this command from the root folder of the cloned POMSimulator repository.

# Route 1: Standard Python Installation with venv and pip

## Install Python 3.10, 3.11, or 3.12

POMSimulator should be installed with Python 3.10, 3.11, or 3.12. Python, 3.12 is a good default choice.

Avoid installing POMSimulator directly into the system's Python environment. Always create a virtual environment.

### Linux

On Ubuntu or Debian:

- sudo apt update
- sudo apt install python3 python3-venv python3-pip python3-dev build-essential


Check the installed Python version:

python3 --version


The version should be 3.10, 3.11, or 3.12.

### macOS

Install Python from:

[https://www.python.org/downloads/](https://www.python.org/downloads/)

Or using Homebrew:

- brew install python


Check the version:

- python3 --version


### Windows

Download Python from:

[https://www.python.org/downloads/windows/](https://www.python.org/downloads/windows/)

Enable `Add python.exe to PATH` during installation.

Check the version:

python --version


## Install Git

Git is required to clone the repository and switch to the `gui\_alpha` branch.

### Linux

sudo apt update

sudo apt install git

git --version


### macOS

brew install git


or

xcode-select --install


Check:

git --version


### Windows

Download Git:

[https://git-scm.com/download/win](https://git-scm.com/download/win)

Check:

git --version


## Clone POMSimulator and Switch to the gui\_alpha Branch

Using the public repository:

git clone https://github.com/petrusen/pomsimulator.git
    cd pomsimulator
    git fetch --all
    git checkout gui\_alpha
    git pull origin gui\_alpha


Verify the active branch:

git branch --show-current


Expected output:

gui\_alpha


## Create and Activate a Virtual Environment

### Linux or macOS

- python3 -m venv .venv
- source .venv/bin/activate
- python -m pip install --upgrade pip setuptools wheel


### Windows PowerShell

- python -m venv .venv
    .\\.venv\\Scripts\\Activate.ps1
- python -m pip install --upgrade pip setuptools wheel


If activation is blocked:

- Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
- .\\.venv\\Scripts\\Activate.ps1


### Windows Command Prompt

- python -m venv .venv
    .venv\\Scripts\\activate.bat
- python -m pip install --upgrade pip setuptools wheel


## Install POMSimulator

Make sure:

* You are inside the cloned repository.
* The active branch is `gui_alpha`.
* The virtual environment is activated.

Install in editable mode:

- pip install -e .

Dependencies include:

* NumPy
* Matplotlib
* NetworkX
* SciPy
* Pandas
* Scikit-learn
* Seaborn
* PyQt5
* ASE

## Verify the Installation

- python --version
- python -m pip --version
- python -c "import pomsimulator; print('POMSimulator import OK')"
- python -c "import PyQt5; print('PyQt5 import OK')"
- python -c "import ase; print('ASE import OK')"


Verify the import location:

python -c "import pomsimulator, pathlib; \\
    print(pathlib.Path(pomsimulator.\_\_file\_\_).resolve())"


## Run the GUI

Linux/macOS:

- cd /path/to/pomsimulator
-    source .venv/bin/activate
-    python GUI.py


Windows PowerShell:

- cd C:\\path\\to\\pomsimulator
- .\\.venv\\Scripts\\Activate.ps1
-    python GUI.py


Windows CMD:

- cd C:\\path\\to\\pomsimulator
- .venv\\Scripts\\activate.bat
-    python GUI.py


## Minimal Command Summary

Linux/macOS:

- git clone https://github.com/petrusen/pomsimulator.git
 -   cd pomsimulator
 -   git checkout gui\_alpha
 -   python3 -m venv .venv
 -   source .venv/bin/activate
 -   python -m pip install --upgrade pip setuptools wheel
 -   pip install -e .
 -   python GUI.py


Windows PowerShell:

- git clone https://github.com/petrusen/pomsimulator.git
-    cd pomsimulator
-    git checkout gui\_alpha
-    python -m venv .venv
-    .\\.venv\\Scripts\\Activate.ps1
-    python -m pip install --upgrade pip setuptools wheel
-    pip install -e .
-    python GUI.py


# Route 2: Conda or Miniconda Installation

## Install Miniconda

Download from:

[https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)

*Make sure to add to path during installation*
Verify installation:

- conda --version


## Clone Repository

- git clone https://github.com/petrusen/pomsimulator.git
-    cd pomsimulator
-    git fetch --all
-    git checkout gui_alpha
-    git pull origin gui_alpha


## Create and Activate Conda Environment

Python 3.11:

- conda create -n pomsimulator_gui python=3.12 pip -y
- conda activate pomsimulator_gui


Alternative versions:

- conda create -n pomsimulator_gui python=3.10 pip -y


or

- conda create -n pomsimulator\_gui python=3.11 pip -y


Upgrade packaging tools:

- python -m pip install --upgrade pip setuptools wheel


## Install POMSimulator

pip install -e .


## Verify Installation

- conda info --envs
-    python --version
-    python -m pip --version
-    python -c "import pomsimulator; print('POMSimulator import OK')"
-    python -c "import PyQt5; print('PyQt5 import OK')"
-    python -c "import ase; print('ASE import OK')"


## Run the GUI

Linux/macOS:

- cd /path/to/pomsimulator
-    conda activate pomsimulator\_gui
-    python GUI.py


Windows:

- cd C:\\path\\to\\pomsimulator
-    conda activate pomsimulator\_gui
-    python GUI.py


## Minimal Conda Summary

- git clone https://github.com/petrusen/pomsimulator.git
-    cd pomsimulator
-    git checkout gui\_alpha
-    conda create -n pomsimulator\_gui python=3.11 pip -y
- conda activate pomsimulator\_gui
-   python -m pip install --upgrade pip setuptools wheel
 -  pip install -e .
 -  python GUI.py



# Final Installation Checklist

The installation is complete when:

* Python 3.10, 3.11, or 3.12 is installed.
* Git is installed.
* The POMSimulator repository is cloned.
* The active branch is `gui_alpha`.
* A virtual environment or Conda environment is activated.
* `pip install -e .` completes successfully.
* `pomsimulator`, `PyQt5`, and `ase` import successfully.
* The GUI opens with `python GUI.py`.

# Troubleshooting \& Error Messages

For any problems or incidents during the installation, don’t hesitate to
contact through email: jbuils@iciq.es

