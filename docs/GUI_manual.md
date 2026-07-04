# POMSimulator GUI — User Manual (Markdown)

> **Version:** 2.0 · **License:** AGPL-3.0 · **Authors:** Enric Petrus, Mireia Segado‑Centellas, Carles Bo · **Developers:** Enric Petrus, Jordi Buils, Diego Garay‑Ruiz.   
> **This manual documents the GUI defined in** [pomsimulator/GUI.py](pomsimulator/GUI.py:1).

---

## Table of Contents

- [1. What is POMSimulator?](#1-what-is-pomsimulator)
- [2. System Requirements](#2-system-requirements)
- [3. Application Layout & UI Basics](#3-application-layout--ui-basics)
  - [3.1 Menus](#31-menus)
  - [3.2 Modes (IPA/HPA)](#32-modes-ipahpa)
  - [3.3 File Dock (Project Browser)](#33-file-dock-project-browser)
  - [3.4 Console & Logs](#34-console--logs)
  - [3.5 Keyboard Shortcuts](#35-keyboard-shortcuts)
- [4. Tabs & Workflows](#4-tabs--workflows)
  - [4.1 Presimulation](#41-presimulation)
    - [4.1.1 Generate Molfiles](#411-generate-molfiles)
    - [4.1.2 Compute Isomorphism](#412-compute-isomorphism)
  - [4.2 Simulation](#42-simulation)
    - [4.2.1 Compute lgkf](#421-compute-lgkf)
  - [4.3 Scaling](#43-scaling)
  - [4.4 Speciation](#44-speciation)
    - [4.4.1 Speciation Diagram](#441-speciation-diagram)
    - [4.4.2 Phase Diagram](#442-phase-diagram)
  - [4.5 Clustering](#45-clustering)
    - [4.5.1 Clustering](#451-clustering)
    - [4.5.2 Model Selection](#452-model-selection)
    - [4.5.3 Boxplot Filtering](#453-boxplot-filtering)
  - [4.6 Plotting](#46-plotting)
    - [4.6.1 Plot Speciation Diagram](#461-plot-speciation-diagram)
    - [4.6.2 Plot Phase Diagram](#462-plot-phase-diagram)
- [5. Q&A / Troubleshooting](#5-qa--troubleshooting)
- [6. Documentation & References](#6-documentation--references)

---

## 1. What is POMSimulator?

**POMSimulator** is a Python framework to predict aqueous speciation and self‑assembly of polyoxometalates (POMs). It builds reaction maps via graph theory, solves the non‑linear speciation model, and outputs (i) a chemical reaction network and (ii) predicted formation constants. Utilities for scaling and plotting (speciation/phase diagrams) are provided. 

This GUI wraps those workflows—reading ADF outputs (or user‑provided data), running simulations, scaling constants, generating speciation/phase diagrams, and analyzing clusters—into one interface. 

---

## 2. System Requirements

**Tested platform / Python:** Python **3.8–3.12** on **Ubuntu 22.04**.

**Minimum (small projects / basic GUI use)**
- CPU: 4 cores  
- RAM: 8 GB  
- Disk: 5 GB free  
- GPU: not required

**Recommended (larger runs / clustering)**
- CPU: 8–16 cores (GUI exposes *cores* for parallel work)   
- RAM: 16–32 GB  
- Disk: 20+ GB free

> *Notes*  
> • Long tasks (e.g., isomorphism, simulation, clustering) run in worker threads and stream logs to the console. More cores improve throughput.   
> • Viewing `.mol` files uses an internal viewer; if an ASE‑based plugin is missing, the GUI will notify you.
---

## 3. Application Layout & UI Basics

The main window shows a **mode toolbar**, the **menu bar**, a **left file dock**, **tabbed workflows** on the right, and a **shared console** at the bottom. 

### 3.1 Menus

**File**
- **Import Configuration File (Ctrl+O)** — Load a `.pomsim` configuration.
- **Export Configuration File... (Ctrl+S)** — Save current parameters to `.pomsim`.
- **Exit (Ctrl+Q)** — Quit the application.

**View**
- **Show Hidden Files** — Toggle visibility of hidden files/folders in the File System dock.

**Tools**
- **Preferences (Ctrl+P)** — Theme (Light/Dark), global text size.
- **Toggle Theme (Ctrl+T)** — Switch Light/Dark.
- **Toggle Mode (Ctrl+M)** — Switch **IPA** ↔ **HPA**.
- **Refresh GUI (Ctrl+F5)** — Reset parameters & indicators.

**Help**
- **User Manual (F1)** — Open the embedded manual.
- **Online Documentation (F2)** — Open ReadTheDocs in your browser.
- **About (Ctrl+A)** — Shows version/authors/license.
- **Key Bindings (Ctrl+F1)** — Full shortcut list.

### 3.2 Modes (IPA/HPA)

Use the toolbar buttons (**IPA**/**HPA**) or **Ctrl+M** to change mode. Switching mode rebuilds mode‑dependent tabs/fields and updates the window title and indicators to prevent mixing IPA/HPA parameters. 

### 3.3 File Dock (Project Browser)

The **File System** dock (left) displays your project tree (works in both source and PyInstaller builds). Double‑click to open:  
- Text (`.txt`, `.py`, `.csv`, `.json`, `.xml`, `.md`, `.log`, `.ini`, `.pomsim`, `.out`) → text viewer (Ctrl+S to save edits).  
- Images (`.png`, `.jpg`, `.jpeg`, `.gif`, `.bmp`, `.tiff`, `.svg`) → image viewer.  
- Molecules (`.mol`) → molecule viewer.  
- Others → open with your OS default app. 

### 3.4 Console & Logs

A bottom **Console Output** panel streams logs from every long‑running task (generation, isomorphism, simulation, scaling, speciation, clustering). Use **Clear** to wipe its contents. Theme/font size propagate here as well. 

### 3.5 Keyboard Shortcuts

- **Ctrl+O/S/Q** — import/export config, exit
- **Ctrl+P/T/M/F5** — preferences, theme, mode, refresh
- **Ctrl++ / Ctrl+-** — increase/decrease text size
- **Ctrl+F1** — key bindings dialog; **F1** — user manual; **F2** — online documentation; **Ctrl+A** — about.

---

## 4. Tabs & Workflows

Each main tab may contain subtabs. Status indicators (✓ / ✗ / ♞ / ∘) appear on subtab headers and on the parent tab. 

### 4.1 Presimulation

Subtabs: **Generate Molfiles**, **Compute Isomorphism**. 

#### 4.1.1 Generate Molfiles

**Purpose** — Create `.mol` files from ADF input/output folders. 

**Inputs**
- **ADF Folder** — directory with ADF outputs. *(Browse)*  
- **Molfile Directory** — destination for generated `.mol`. *(Browse)* 

**Run / Stop**
- **▶ Run Generate Molfiles** — launches the generator in a background worker; logs stream to the console.  
- **⏹ Stop** — requests graceful termination. 

**Validation**
- If either folder is missing you’ll see: **“Please select both ADF and Molfile directories.”** 

#### 4.1.2 Compute Isomorphism

**Purpose** — Build an isomorphism matrix across generated molecules. 

**Inputs**
- **System** — POM system name (e.g., `W`, `Mo`, …).  
- **Molfile Directory** — folder with `.mol` files. *(Browse)*  
- **Output Directory** — destination for isomorphism data. *(Browse)*  
- **Isomorphism cores** — number of CPU cores to use. 

**Run / Stop**
- **▶ Run Compute Isomorphism** — executes in a worker thread; progress logs stream to the console.  
- **⏹ Stop** — cancels gracefully. 

**Validation**
- Missing any required path or the system name triggers the same **ADF/Molfile** warning. 

---

### 4.2 Simulation

Subtab: **Compute lgkf** (formation constants). 

#### 4.2.1 Compute lgkf

**Purpose** — Compute formation constants for the chosen system and conditions. Mode‑aware (**IPA** or **HPA**). 

**Required preparation**
- **POM System**, **ADF Folder**, **Molfile Directory**, **Output Path**.  
- **Reference compound(s)** — at least one (two in HPA). 

**General settings** (typical fields)
- **Cores**, **Batch size**, **Sample %**, **Sample type** (`random`/`all`).
- Chemical settings such as **pH window (min/max/step)**, **Temperature (K)**, **Ionic strength**, **Energy threshold**, and **Proton difference threshold**.

**Mode‑specific inputs**
- **IPA** — **C₀** (total concentration, M) and **Reference compound** label.  
- **HPA** — **C_M**, **C_X** (reagent totals), **M** and **X** reference compound labels. 

**Run / Stop**
- **▶ Run Simulation** — starts the appropriate backend in a worker; logs appear in the console.  
- **⏹ Stop** — cancels gracefully. 

---

### 4.3 Scaling

**Purpose** — Linearly scale predicted formation constants to match selected experimental sets before plotting speciation/phase diagrams. 

**Inputs**
- **POM System** (text) and **Output Path** (directory). *(Browse)*  
- **Scaling Mode** — options populated from the internal database.  
- **Experimental Set** — choose a pre‑encoded set; the GUI displays values for reference. 

**Run / Stop**
- **▶ Run Scaling** — executes in a worker; results are saved under the chosen output path.  
- **⏹ Stop** — cancels gracefully. 

---

### 4.4 Speciation

Subtabs: **Speciation Diagram**, **Phase Diagram** (both are mode‑aware). 

#### 4.4.1 Speciation Diagram

**Purpose** — Compute and plot species distribution vs pH from scaled formation constants.

**Prerequisites**
- Requires `logkf_<SYSTEM>.csv` and `scaling_params_<SYSTEM>.pomsim` under the selected **Output Path** (produced by Simulation and Scaling).

**Inputs**
- **System** and **Output Path**. *(Browse)*
- **Species to Plot** — Labels File and Selected Labels (or `all`).
- **pH Range** — min/max/step.
- **Mode‑specific**:
  - **IPA** — Initial Concentration (**C**, mol/L).
  - **HPA** — Initial Heteroatom Concentration (**C_X**, mol/L) and Initial Metal Concentration (**C_M**, mol/L).
- **Operation** — **Cores**, **Batch size**.

**Run / Stop**
- **▶ Generate Speciation Diagram** — saves `Array_<SYSTEM>.npz`, `speciation_params_<SYSTEM>.txt`, and figures under the **Output Path**; logs to console.
- **⏹ Stop** — cancels gracefully.

#### 4.4.2 Phase Diagram

**Purpose** — Generate phase diagrams (e.g., concentration‑/pH‑based maps).

**Inputs (mode‑specific)**
- **System** and **Output Path**. *(Browse)*
- **Phase Diagram Directory Name** — target folder name for the outputs.
- **Model Subset File** — optional file to restrict simulated models.
- **Species to Plot** — Labels File and Selected Labels (or `all`).
- **pH Range** — min/max/step.
- **IPA**: concentration grid — Min/Max Concentration (mol/L) and Number of points.
- **HPA**: ratio grid — Min/Max Metal/Heteroatom ratio and Number of points; Initial Heteroatom Concentration (**C_X**, mol/L).
- **Operation** — **Cores**, **Batch size**.
- **Outputs** — phase directory with `array_XX.npz` files and `npz_info_<SYSTEM>.dat` under the **Output Path**.

**Run / Stop**
- **▶ Generate Phase Diagram** — runs in a worker and saves figures under the output path.
- **⏹ Stop** — cancels gracefully.

---

### 4.5 Clustering

Subtabs: **Clustering**, **Model Selection**, **Boxplot Filtering**. 

#### 4.5.1 Clustering

**Purpose** — Compute feature‑based clusters and figures for selected species. 

**Inputs**
- **System** (text), **Output Path** (directory). *(Browse)*  
- **Cluster Directory / NPZ file** — where to read/write cluster data.  
- **Features file** — numeric features per species; **Feature Selection** from a built‑in feature list.  
- **n_clusters** — number of clusters; **normalize features** toggle.  
- **Labels / Color Dictionary** — select species labels and an optional color map for consistent plots. 

**Run / Stop**
- **▶ Run Clustering** — launches clustering; figures and data are saved in the cluster directory.  
- **⏹ Stop** — cancels gracefully. 

#### 4.5.2 Model Selection

**Purpose** — Post‑hoc selection of cluster(s) and features for downstream analysis/plots.  
- **Cluster Selection** — check boxes `0..n_clusters‑1`.  
- **Feature Selection** — same feature list as in **Clustering** for quick refinement. 

#### 4.5.3 Boxplot Filtering

**Purpose** — Choose metrics to visualize as boxplots and filter groups accordingly.  
- Configure **boxplot list** and **color dictionary**; saved under **Visualization** in the configuration. 

### 4.6 Plotting

Subtabs: **Plot Speciation Diagram**, **Plot Phase Diagram**.

#### 4.6.1 Plot Speciation Diagram

**Purpose** — Plot speciation diagrams from previously computed NPZ arrays.

**Inputs**
- **System** — POM system name.
- **Output Path** — directory containing results.
- **NPZ File** — path to `Array_<SYSTEM>.npz` or a compatible NPZ file.
- **Metal Index (m_idx)** — `0` for IPA; `0` or `1` for HPA.
- **Species to Plot** — Labels File and Selected Labels (or `all`).
- **Color Dictionary** — optional predefined color map.

**Run / Stop**
- **▶ Generate Speciation Plot** — creates `Speciation_Diagram_<SYSTEM>.png` under the output path; logs stream to the console.
- **⏹ Stop** — requests graceful termination.

#### 4.6.2 Plot Phase Diagram

**Purpose** — Plot phase diagrams from phase calculation results.

**Inputs**
- **System** — POM system name.
- **Output Path** — directory containing results.
- **Phase Diagram Directory Name** — folder containing `array_XX.npz` files and `npz_info_<SYSTEM>.dat`.
- **Species to Plot** — Labels File and Selected Labels (or `all`).
- **Color Dictionary** — optional predefined color map.

**Run / Stop**
- **▶ Generate Phase Plot** — writes figures under the output path; logs stream to the console.
- **⏹ Stop** — requests graceful termination.

---

## 5. Q&A / Troubleshooting

**Q1. The GUI warns: “Please select both ADF and Molfile directories.”**  
**A.** Provide valid folders in **Presimulation → Generate Molfiles** and/or **Compute Isomorphism**. 

**Q2. Molecule viewer fails to open a `.mol` file.**  
**A.** Ensure the file is valid and that the molecule viewer dependencies are available; the GUI will display a helpful message if an importer is missing. 

**Q3. After switching IPA/HPA, my fields reset.**  
**A.** That’s expected—mode toggling rebuilds mode‑dependent widgets to avoid mixing IPA/HPA parameters. 

**Q4. Where can I see progress and errors?**  
**A.** The **Console Output** panel shows timestamped logs for all operations; **Clear** wipes the panel. 

---

## 6. Documentation & References

### Official docs & resources
- **ReadTheDocs (latest):** https://pomsimulator.readthedocs.io/en/latest/ (also available via **Help → Documentation**). 

### How to cite POMSimulator
When publishing results generated with POMSimulator, please cite: 

- Petrus, E.; Segado, M.; Bo, C. *Chem. Sci.* **2020**, 11, 8448–8456. DOI: 10.1039/D0SC03530K  
- Petrus, E.; Buils, J.; Garay‑Ruiz, D.; Segado‑Centellas, M.; Bo, C. *J. Comput. Chem.* **2024**, 45, 2242–2250. DOI: 10.1002/jcc.27389

**Selected peer‑reviewed articles featuring POMSimulator:** 
- Petrus, E.; Bo, C. *J. Phys. Chem. A.* **2021**, 125, 5212–5219.  
- Petrus, E.; Segado‑Centellas, M.; Bo, C. *Inorg. Chem.* **2022**, 61, 13708–13718.  
- Petrus, E.; Garay‑Ruiz, D.; Reiher, M.; Bo, C. *J. Am. Chem. Soc.*, **2023**, 145, 18920–18930.  
- Garay‑Ruiz, D.; Petrus, E.; Bo, C. *Rev. Soc. Catalana de Química* **2023**, 22, 23–38.  
- Buils, J.; Garay‑Ruiz, D.; Segado‑Centellas, M.; Petrus, E.; Bo, C. *Chem. Sci.* **2024**, 15, 14218–14227.  
- Buils, J.; Garay‑Ruiz, D.; Petrus, E.; Segado‑Centellas, M.; Bo, C. *Digital Discovery* **2025**, 4, 970–978.

---

### Quick Start (suggested order)

1) **Presimulation → Generate Molfiles**: pick **ADF Folder** and **Molfile Directory** → **Run**.   
2) **Presimulation → Compute Isomorphism**: set **System**, **Molfile Directory**, **Output Directory**, **Cores** → **Run**.   
3) **Simulation → Compute lgkf**: provide **Preparation** paths + **Reference compound(s)**; set pH/I/T; choose solver → **Run**.   
4) **Scaling**: choose **Scaling Mode** and **Experimental Set**, ensure **Output Path** → **Run Scaling**.   
5) **Speciation**: select labels/colors; **Speciation Diagram** or **Phase Diagram** → **Generate**.   
6) **Clustering** (optional): pick features, set **n_clusters**, labels/colors → **Run Clustering**, then refine via **Model Selection** and **Boxplot Filtering**. 

---

*This manual summarizes the behavior of the GUI implementation.*
