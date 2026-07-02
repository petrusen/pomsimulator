# POMSimulator User Manual

**Version 1.0**  
**Date: July 2026**

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [System Requirements](#2-system-requirements)
3. [Installation](#3-installation)
4. [Getting Started](#4-getting-started)
5. [Interface Overview](#5-interface-overview)
6. [File System Navigation](#6-file-system-navigation)
7. [Preview System](#7-preview-system)
8. [Menu System](#8-menu-system)
9. [Simulation Workflow](#9-simulation-workflow)
10. [Presimulation Tab](#10-presimulation-tab)
11. [Simulation Tab](#11-simulation-tab)
12. [Scaling Tab](#12-scaling-tab)
13. [Speciation Tab](#13-speciation-tab)
14. [Clustering Tab](#14-clustering-tab)
15. [Plotting Tab](#15-plotting-tab)
16. [Configuration Management](#16-configuration-management)
17. [Console Output](#17-console-output)
18. [Keyboard Shortcuts](#18-keyboard-shortcuts)
19. [Troubleshooting](#19-troubleshooting)
20. [Frequently Asked Questions](#20-frequently-asked-questions)

---

## 1. Introduction

POMSimulator is a comprehensive scientific software application designed for the simulation and analysis of Polyoxometalate (POM) systems. The application provides a complete workflow for:

- **Molecular Structure Generation**: Creating molecular structure files from quantum chemistry calculations
- **Isomorphism Analysis**: Computing structural similarities between molecules
- **Formation Constant Calculations**: Determining thermodynamic stability parameters
- **Scaling Analysis**: Calibrating computational results with experimental data
- **Speciation Modeling**: Predicting species distribution in solution
- **Clustering Analysis**: Grouping similar molecular structures
- **Data Visualization**: Creating publication-quality plots and diagrams

The software is built using PyQt5 and provides an intuitive graphical user interface that guides users through the complete simulation workflow from initial data preparation to final results visualization.

### 1.1 Target Audience

This manual is intended for:
- Computational chemists working with polyoxometalate systems
- Researchers studying solution chemistry and speciation
- Graduate students learning molecular simulation techniques
- Scientists requiring thermodynamic modeling capabilities

### 1.2 Prerequisites

Users should have:
- Basic understanding of quantum chemistry concepts
- Familiarity with molecular structure formats (MOL files, ADF output files)
- Knowledge of thermodynamic principles
- Experience with scientific data analysis

---

## 2. System Requirements

### 2.1 Minimum Requirements

- **Operating System**: Windows 10/11, macOS 10.14+, or Linux (Ubuntu 18.04+)
- **Memory**: 4 GB RAM minimum, 8 GB recommended
- **Storage**: 2 GB available disk space
- **Python**: Python 3.10, 3.11 or 3.12

### 2.2 Required Dependencies

- PyQt5 5.12+
- NumPy 1.18+
- Pandas 1.0+
- Matplotlib 3.1+
- SciPy 1.4+
- ASE (Atomic Simulation Environment) 3.19+
- NetworkX 2.4+
- Scikit-learn 0.22+

### 2.3 Optional Dependencies

- Jupyter Notebook (for advanced analysis)
- OpenGL drivers (for 3D molecular visualization)

---

## 3. Installation

### 3.1 Installing from Source

1. **Download the source code** from the repository
2. **Navigate to the project directory**:
   ```bash
   cd pomsimulator
   ```
3. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
4. **Run the application**:
   ```bash
   python GUI.py
   ```

### 3.2 Environment Setup

For optimal performance, create a dedicated Python environment:

```bash
conda create -n pomsimulator python=3.8
conda activate pomsimulator
pip install -r requirements.txt
```

### 3.3 Verifying Installation

Launch the application and verify that:
- The main window opens without errors
- All tabs are accessible
- The file browser displays the project directory
- Sample data files are visible in the inputs folder

---

## 4. Getting Started

### 4.1 First Launch

When you first launch POMSimulator:

1. **Main Window**: The application opens with a large main window containing multiple tabs
2. **File Browser**: The left panel shows a file system navigator
3. **Console**: The bottom panel displays system messages and output
4. **Tabs**: Six main tabs organize different workflow stages

### 4.2 Quick Start Workflow

For new users, follow this basic workflow:

1. **Load Sample Data**: Navigate to the `inputs` folder and explore sample configurations
2. **Import Configuration**: Use `File > Import Configuration File` to load a sample `.pomsim` file
3. **Review Parameters**: Check the loaded parameters in each tab
4. **Run Simulation**: Start with the Presimulation tab and work through each stage
5. **View Results**: Check the `outputs` folder for generated files and plots

### 4.3 Understanding the Interface

The interface is organized into logical sections:
- **Left Panel**: File navigation and preview area
- **Center Panel**: Main workflow tabs
- **Bottom Panel**: Console output and messages
- **Top Panel**: Menu bar and toolbar

---

## 5. Interface Overview

### 5.1 Main Window Layout

> **Figure 1.** Main Window Layout

The main window uses a three-panel layout:

**Left Panel (File System Dock)**:
- File browser for project navigation
- Preview area for images and molecules
- Resizable splitter between browser and preview

**Center Panel (Tab Widget)**:
- Six main workflow tabs
- Each tab contains specific simulation tools
- Progress indicators show completion status

**Bottom Panel (Console)**:
- Real-time output from operations
- Error messages and warnings
- Execution logs and progress updates

### 5.2 Window Controls

**Resizing**: All panels can be resized by dragging the splitter bars
**Docking**: The file system panel can be moved to the right side
**Themes**: Switch between light and dark themes via the Tools menu
**Scaling**: Interface elements automatically scale with system DPI settings

### 5.3 Status Indicators

The application provides visual feedback through:
- **Tab Colors**: Indicate completion status (idle, running, completed, error)
- **Button States**: Show available actions (enabled/disabled)
- **Progress Bars**: Display operation progress where applicable
- **Console Messages**: Provide detailed status information

---

## 6. File System Navigation

### 6.1 File Browser Panel

> **Figure 2.** File Browser Panel

The file browser panel provides access to all project files:

**Navigation Features**:
- **Tree View**: Hierarchical display of folders and files
- **File Types**: Icons indicate different file types
- **Sorting**: Files are sorted alphabetically
- **Expansion**: Click folder icons to expand/collapse directories

**Supported File Types**:
- **Configuration Files** (`.pomsim`): Project settings and parameters
- **Output Files** (`.out`): Quantum chemistry calculation results
- **Molecule Files** (`.mol`): 3D molecular structure data
- **Image Files** (`.png`, `.jpg`, `.svg`): Plots and diagrams
- **Data Files** (`.csv`, `.txt`, `.npz`): Numerical results and arrays
- **Text Files** (`.py`, `.md`, `.log`): Scripts and documentation

### 6.2 File Operations

**Double-Click Actions**:
- **Text Files**: Open in built-in text editor
- **Image Files**: Display in preview panel with zoom controls
- **Molecule Files**: Show 3D structure in molecule viewer
- **Other Files**: Open with system default application

**Context Menu** (Right-click):
- Copy file path
- Open containing folder
- Refresh directory view

### 6.3 Directory Structure

The typical project structure includes:

```
pomsimulator/
├── inputs/                 # Input data and configurations
│   ├── *.pomsim           # Configuration files
│   ├── *_Set/             # Quantum chemistry output files
│   └── *_molfiles/        # Molecular structure files
├── outputs/               # Generated results
│   ├── *_data/           # Simulation results by system
│   └── *.png             # Generated plots
├── docs/                  # Documentation
└── utilities/             # Helper scripts
```

---

## 7. Preview System

### 7.1 Preview Manager

> **Figure 3.** Preview System with Image Viewer

The preview system allows you to view files without opening separate windows:

**Preview Types**:
- **Image Preview**: For plots, diagrams, and figures
- **Molecule Preview**: For 3D molecular structures
- **No Preview**: When no compatible file is selected

### 7.2 Image Viewer

**Features**:
- **Zoom Controls**: Zoom in/out with 25% increments
- **Pan Support**: Drag to move around zoomed images
- **Navigation**: Browse all images in the same folder
- **Status Display**: Shows filename, image count, and zoom level

**Controls**:
- **Previous/Next**: Navigate between images in folder
- **Zoom In (+)**: Increase image size by 25%
- **Zoom Out (-)**: Decrease image size by 25%
- **Reset Zoom**: Return to 100% size
- **Close**: Hide the preview panel

**Mouse Interactions**:
- **Left-click and drag**: Pan around zoomed images
- **Scroll wheel**: Zoom in/out (if supported)

### 7.3 Molecule Viewer

> **Figure 4.** 3D Molecule Viewer

**Features**:
- **3D Visualization**: Interactive 3D molecular structures
- **Rotation**: Click and drag to rotate molecules
- **Element Colors**: Standard CPK coloring scheme
- **Bond Display**: Shows chemical bonds between atoms

**Controls**:
- **Rotation**: Left-click and drag to rotate view
- **Zoom**: Mouse wheel or zoom controls
- **Reset View**: Return to default orientation
- **Close**: Hide the molecule viewer

### 7.4 Preview Panel Management

**Automatic Sizing**: The preview panel automatically adjusts its size when content is loaded
**Collapsible**: Close previews to maximize file browser space
**Persistent**: Preview settings are maintained during the session

---

## 8. Menu System

### 8.1 File Menu

> **Figure 5.** File Menu Options

**Import Configuration File (Ctrl+O)**:
- Load simulation parameters from `.pomsim` files
- Populates all tabs with saved settings
- Validates parameter compatibility

**Export Configuration File (Ctrl+S)**:
- Save current parameters to `.pomsim` format
- Preserves all tab settings
- Creates portable configuration files

**Exit (Ctrl+Q)**:
- Close the application
- Prompts to save unsaved work

### 8.2 Tools Menu

**Preferences (Ctrl+P)**:
- Open application settings dialog
- Configure default parameters
- Adjust interface preferences

**Toggle Theme (Ctrl+T)**:
- Switch between light and dark themes
- Applies immediately to entire interface
- Preference is saved between sessions

**Toggle Mode (Ctrl+M)**:
- Switch between simulation and visualization modes
- Changes available functionality
- Affects tab visibility and options

**Refresh GUI (Ctrl+F5)**:
- Reset interface to default state
- Clear all parameters
- Useful for starting fresh

### 8.3 Help Menu

**Documentation (F2)**:
- Open online documentation
- Access user guides and tutorials
- View API reference

**About (Ctrl+A)**:
- Display application information
- Show version number and credits
- License information

**Key Bindings (Ctrl+F1)**:
- Display keyboard shortcuts reference
- Organized by category
- Searchable shortcut list

---

## 9. Simulation Workflow

### 9.1 Overview

POMSimulator follows a structured workflow with six main stages:

1. **Presimulation**: Prepare molecular data and compute structural properties
2. **Simulation**: Calculate formation constants and thermodynamic parameters
3. **Scaling**: Calibrate computational results with experimental data
4. **Speciation**: Model species distribution in solution
5. **Clustering**: Analyze and group molecular structures
6. **Plotting**: Generate publication-quality visualizations

### 9.2 Workflow Dependencies

Each stage depends on the completion of previous stages:

```
Presimulation → Simulation → Scaling → Speciation
                                    ↓
                            Clustering → Plotting
```

**Sequential Processing**: Complete stages in order for best results
**Parallel Branches**: Clustering and Speciation can run independently after Scaling
**Iterative Refinement**: Return to earlier stages to refine parameters

### 9.3 Data Flow

**Input Data**:
- Quantum chemistry output files (`.out`)
- Experimental formation constants
- Molecular structure files (`.mol`)

**Intermediate Data**:
- Computed formation constants
- Scaling parameters
- Feature matrices

**Output Data**:
- Speciation diagrams
- Phase diagrams
- Clustering results
- Statistical analyses

### 9.4 Progress Tracking

**Tab Indicators**: Visual indicators show completion status
**Console Output**: Detailed progress messages
**File Generation**: Output files confirm successful completion
**Error Handling**: Clear error messages guide troubleshooting

---

## 10. Presimulation Tab

### 10.1 Overview

> **Figure 6.** Presimulation Tab Interface

The Presimulation tab prepares molecular data for subsequent analysis. It contains two main subtabs:

1. **Generate Molfiles**: Convert quantum chemistry outputs to molecular structure files
2. **Compute Isomorphism**: Analyze structural similarities between molecules

### 10.2 Generate Molfiles Subtab

**Purpose**: Convert ADF (Amsterdam Density Functional) output files to MOL format for molecular visualization and analysis.

**Parameters**:

**ADF Folder**:
- **Description**: Directory containing quantum chemistry output files
- **Format**: `.out` files from ADF calculations
- **Selection**: Use Browse button to select folder
- **Validation**: Folder must contain valid `.out` files

**Molfile Directory**:
- **Description**: Output directory for generated MOL files
- **Format**: Standard MOL file format
- **Creation**: Directory is created if it doesn't exist
- **Organization**: Files are named based on input file names

**Workflow**:
1. **Select Input Folder**: Browse to folder containing `.out` files
2. **Choose Output Location**: Specify where MOL files will be saved
3. **Run Generation**: Click "▶ Run Generate Molfiles"
4. **Monitor Progress**: Watch console for processing updates
5. **Verify Results**: Check output folder for generated MOL files

**Output Files**:
- One `.mol` file per input `.out` file
- Files contain 3D atomic coordinates
- Compatible with molecular visualization software

### 10.3 Compute Isomorphism Subtab

**Purpose**: Analyze structural relationships between molecules to identify similar species and build reaction networks.

**Parameters**:

**Input Directory**:
- **Description**: Folder containing MOL files for analysis
- **Requirements**: MOL files from Generate Molfiles step
- **Validation**: Must contain valid molecular structures

**Analysis Options**:
- **Similarity Threshold**: Minimum similarity for isomorphism detection
- **Bond Tolerance**: Allowed variation in bond lengths
- **Angle Tolerance**: Allowed variation in bond angles

**Workflow**:
1. **Select MOL Directory**: Choose folder with molecular structure files
2. **Configure Parameters**: Set similarity thresholds and tolerances
3. **Run Analysis**: Click "▶ Run Compute Isomorphism"
4. **Review Results**: Examine isomorphism matrix and groupings
5. **Export Data**: Save results for use in subsequent steps

**Output Files**:
- Isomorphism matrix showing structural similarities
- Grouping information for similar molecules
- Network topology data

### 10.4 Best Practices

**File Organization**:
- Keep input and output folders well-organized
- Use descriptive folder names
- Maintain consistent naming conventions

**Parameter Selection**:
- Start with default tolerance values
- Adjust based on chemical system requirements
- Document parameter choices for reproducibility

**Quality Control**:
- Verify MOL file generation completeness
- Check isomorphism results for chemical reasonableness
- Review console output for warnings or errors

---

## 11. Simulation Tab

### 11.1 Overview

> **Figure 7.** Simulation Tab Interface

The Simulation tab performs the core thermodynamic calculations, computing formation constants for polyoxometalate species based on quantum chemistry data.

### 11.2 System Configuration

**POM System Selection**:
- **IPA (Isopolyanions)**: Systems containing only one type of metal (e.g., tungsten, molybdenum)
- **HPA (Heteropolyanions)**: Systems with multiple metal types (e.g., phosphomolybdates)

**System Parameters**:

**Metal Type**:
- **Options**: W (Tungsten), Mo (Molybdenum), As (Arsenic), C (Carbon), P (Phosphorus)
- **Selection**: Choose based on your chemical system
- **Impact**: Determines calculation methods and reference data

**pH Range**:
- **Minimum pH**: Lower bound for calculations (typically 0-2)
- **Maximum pH**: Upper bound for calculations (typically 12-14)
- **Resolution**: Number of pH points for calculations

**Ionic Strength**:
- **Value**: Solution ionic strength in mol/L
- **Range**: Typically 0.1 to 3.0 M
- **Impact**: Affects activity coefficients and equilibrium constants

### 11.3 Computational Parameters

**Energy Calculation Settings**:

**Temperature**:
- **Value**: Calculation temperature in Kelvin
- **Default**: 298.15 K (25°C)
- **Range**: 273-373 K for aqueous systems

**Pressure**:
- **Value**: System pressure in atm
- **Default**: 1.0 atm
- **Impact**: Minor effect for condensed phase systems

**Solvation Model**:
- **COSMO**: Conductor-like Screening Model
- **PCM**: Polarizable Continuum Model
- **Selection**: Based on quantum chemistry method used

### 11.4 Reference Data

**Experimental Constants**:
- **Source**: Literature formation constants
- **Format**: Log K values at specified conditions
- **Validation**: Automatic consistency checking
- **Updates**: Can be modified for different experimental datasets

**Standard States**:
- **Aqueous Species**: 1 M standard state
- **Solid Phases**: Pure solid standard state
- **Gas Phase**: 1 atm standard state

### 11.5 Calculation Workflow

1. **Configure System**: Select POM type and metal composition
2. **Set Conditions**: Define pH range, ionic strength, and temperature
3. **Load Structures**: Import molecular geometries from Presimulation
4. **Run Calculations**: Execute formation constant calculations
5. **Validate Results**: Compare with experimental data where available
6. **Export Data**: Save results for subsequent analysis steps

**Progress Monitoring**:
- Real-time console updates
- Progress bar for long calculations
- Intermediate result validation
- Error detection and reporting

### 11.6 Output Files

**Formation Constants** (`logkf_*.csv`):
- Computed log K values for each species
- Temperature and ionic strength corrections
- Uncertainty estimates where available

**Reaction Network** (`CRN_*.txt`):
- Complete set of formation reactions
- Stoichiometric coefficients
- Reaction pathways and mechanisms

**Simulation Parameters** (`simulation_parameters_*.txt`):
- Complete record of calculation settings
- Reproducibility information
- Version and timestamp data

---

## 12. Scaling Tab

### 12.1 Overview

> **Figure 8.** Scaling Tab Interface

The Scaling tab calibrates computational results with experimental data to improve prediction accuracy. This critical step adjusts theoretical formation constants to match known experimental values.

### 12.2 Scaling Methods

**Linear Scaling**:
- **Method**: Simple linear regression between computed and experimental values
- **Equation**: K_scaled = a × K_computed + b
- **Use Case**: When systematic errors are approximately linear

**Logarithmic Scaling**:
- **Method**: Scaling applied to logarithmic values
- **Equation**: log K_scaled = a × log K_computed + b
- **Use Case**: When errors are proportional to magnitude

**Multi-parameter Scaling**:
- **Method**: Scaling based on molecular descriptors
- **Variables**: Size, charge, coordination number
- **Use Case**: Complex systems with multiple error sources

### 12.3 Reference Dataset Selection

**Experimental Data Sources**:
- **Literature Values**: Peer-reviewed formation constants
- **Database Entries**: Standardized thermodynamic databases
- **Custom Data**: User-provided experimental measurements

**Data Quality Criteria**:
- **Ionic Strength**: Consistent experimental conditions
- **Temperature**: Standard temperature (25°C) preferred
- **pH Range**: Overlapping with computational predictions
- **Uncertainty**: Known error bounds for weighting

### 12.4 Scaling Parameters

**Regression Settings**:

**Fitting Method**:
- **Least Squares**: Standard linear regression
- **Weighted Least Squares**: Account for experimental uncertainties
- **Robust Regression**: Minimize outlier influence

**Cross-Validation**:
- **K-Fold**: Divide data into training and validation sets
- **Leave-One-Out**: Use each point as validation once
- **Bootstrap**: Random sampling with replacement

**Quality Metrics**:
- **R²**: Coefficient of determination
- **RMSE**: Root mean square error
- **MAE**: Mean absolute error

### 12.5 Scaling Workflow

1. **Load Computed Data**: Import formation constants from Simulation tab
2. **Select Reference Data**: Choose experimental dataset for calibration
3. **Configure Method**: Select scaling approach and parameters
4. **Perform Scaling**: Execute calibration calculations
5. **Validate Results**: Assess scaling quality and statistics
6. **Apply Scaling**: Generate scaled formation constants
7. **Export Parameters**: Save scaling coefficients for future use

### 12.6 Quality Assessment

**Statistical Analysis**:
- **Correlation Plots**: Computed vs. experimental values
- **Residual Analysis**: Error distribution and patterns
- **Outlier Detection**: Identification of problematic data points

**Chemical Validation**:
- **Trend Analysis**: Consistency with chemical intuition
- **Systematic Errors**: Detection of method-specific biases
- **Extrapolation Limits**: Valid range for scaled predictions

### 12.7 Output Files

**Scaling Parameters** (`scaling_params_*.pomsim`):
- Calibration coefficients and equations
- Statistical quality metrics
- Validation results and diagnostics

**Scaled Constants** (`scaled_logkf_*.csv`):
- Calibrated formation constants
- Uncertainty propagation
- Confidence intervals

**Regression Analysis** (`regression_output.csv`):
- Detailed statistical analysis
- Residual plots and diagnostics
- Cross-validation results

---

## 13. Speciation Tab

### 13.1 Overview

> **Figure 9.** Speciation Tab Interface

The Speciation tab models the distribution of chemical species in solution as a function of pH, concentration, and other variables. It generates speciation diagrams and phase diagrams for visualization and analysis.

### 13.2 Speciation Calculations

**Species Distribution**:
- **Alpha Diagrams**: Fraction of each species vs. pH
- **Concentration Profiles**: Absolute concentrations vs. conditions
- **Predominance Regions**: Dominant species identification

**Calculation Methods**:
- **Newton-Raphson**: Iterative solution of equilibrium equations
- **Matrix Methods**: Linear algebra approach for large systems
- **Optimization**: Minimization of Gibbs free energy

### 13.3 Diagram Types

**Speciation Diagrams**:
- **X-axis**: pH (typically 0-14)
- **Y-axis**: Species fraction (0-1)
- **Curves**: One curve per chemical species
- **Colors**: Automatic assignment from database

**Phase Diagrams**:
- **2D Plots**: Two variables (e.g., pH vs. concentration)
- **Contour Lines**: Constant species fraction
- **Regions**: Areas of species predominance
- **Boundaries**: Phase transition lines

### 13.4 Parameters

**Solution Conditions**:

**Total Concentration**:
- **Value**: Total metal concentration in mol/L
- **Range**: Typically 10⁻⁶ to 10⁻¹ M
- **Impact**: Affects species distribution patterns

**pH Range**:
- **Minimum**: Lower pH limit (typically 0-2)
- **Maximum**: Upper pH limit (typically 12-14)
- **Resolution**: Number of calculation points

**Ionic Strength**:
- **Value**: Background electrolyte concentration
- **Effect**: Activity coefficient corrections
- **Consistency**: Should match Simulation tab settings

**Temperature**:
- **Value**: System temperature in Kelvin
- **Default**: 298.15 K (25°C)
- **Range**: Limited by thermodynamic data validity

### 13.5 Advanced Options

**Activity Coefficients**:
- **Davies Equation**: Extended Debye-Hückel model
- **Pitzer Model**: Comprehensive ion interaction approach
- **Ideal Solution**: No activity corrections (γ = 1)

**Solid Phases**:
- **Precipitation**: Include solid formation equilibria
- **Solubility**: Automatic solubility limit detection
- **Metastability**: Consider kinetic barriers

**Complexation**:
- **Ligand Binding**: Additional complexing agents
- **Competition**: Multiple ligand systems
- **Side Reactions**: Parallel equilibria

### 13.6 Workflow

1. **Load Formation Constants**: Import scaled data from previous steps
2. **Set Conditions**: Define pH range, concentration, and ionic strength
3. **Configure Options**: Select activity models and additional equilibria
4. **Generate Diagrams**: Calculate and plot speciation distributions
5. **Analyze Results**: Interpret species behavior and transitions
6. **Export Plots**: Save diagrams in various formats
7. **Save Data**: Export numerical results for further analysis

### 13.7 Interpretation Guidelines

**Speciation Diagrams**:
- **Dominant Species**: Highest fraction at given pH
- **Transition Points**: pH where species fractions are equal
- **Buffer Regions**: Areas of gradual species change
- **Sharp Transitions**: Rapid species interconversion

**Phase Diagrams**:
- **Stability Fields**: Regions where species predominate
- **Phase Boundaries**: Lines of equal stability
- **Triple Points**: Intersection of three phase boundaries
- **Critical Points**: Special thermodynamic conditions

### 13.8 Output Files

**Speciation Diagrams** (`Speciation_Diagram_*.png`):
- High-resolution publication-quality plots
- Customizable colors and formatting
- Multiple export formats (PNG, SVG, PDF)

**Numerical Data** (`speciation_data_*.csv`):
- Species fractions vs. pH
- Concentration profiles
- Thermodynamic parameters

**Parameters File** (`speciation_params_*.txt`):
- Complete calculation settings
- Reproducibility information
- Quality metrics and validation

---

## 14. Clustering Tab

### 14.1 Overview

> **Figure 10.** Clustering Tab Interface

The Clustering tab performs unsupervised machine learning analysis to group similar molecular structures and identify patterns in the chemical space. This analysis helps understand structure-property relationships and guide experimental design.

### 14.2 Clustering Methods

**K-Means Clustering**:
- **Algorithm**: Partition data into k clusters
- **Optimization**: Minimize within-cluster variance
- **Parameters**: Number of clusters (k)
- **Use Case**: Well-separated, spherical clusters

**Hierarchical Clustering**:
- **Algorithm**: Build cluster tree (dendrogram)
- **Linkage**: Single, complete, or average linkage
- **Distance**: Various distance metrics available
- **Use Case**: Nested cluster structures

**DBSCAN**:
- **Algorithm**: Density-based clustering
- **Parameters**: Epsilon (neighborhood size), minimum points
- **Advantages**: Handles noise and arbitrary cluster shapes
- **Use Case**: Irregular cluster boundaries

### 14.3 Feature Selection

**Molecular Descriptors**:
- **Geometric**: Bond lengths, angles, volumes
- **Electronic**: Charges, dipole moments, polarizabilities
- **Topological**: Connectivity indices, graph properties
- **Thermodynamic**: Formation energies, stability measures

**Feature Engineering**:
- **Normalization**: Scale features to comparable ranges
- **Selection**: Choose most informative descriptors
- **Dimensionality Reduction**: PCA, t-SNE for visualization
- **Correlation Analysis**: Remove redundant features

### 14.4 Clustering Parameters

**Algorithm Settings**:

**Number of Clusters**:
- **Determination**: Elbow method, silhouette analysis
- **Range**: Typically 2-20 clusters
- **Validation**: Cross-validation and stability analysis

**Distance Metrics**:
- **Euclidean**: Standard geometric distance
- **Manhattan**: Sum of absolute differences
- **Cosine**: Angular similarity measure
- **Custom**: Chemical-specific distance functions

**Preprocessing Options**:
- **Standardization**: Zero mean, unit variance
- **Normalization**: Scale to [0,1] range
- **Outlier Removal**: Detect and exclude anomalous structures

### 14.5 Cluster Analysis

**Quality Metrics**:
- **Silhouette Score**: Cluster separation quality
- **Inertia**: Within-cluster sum of squares
- **Calinski-Harabasz Index**: Cluster validity measure
- **Davies-Bouldin Index**: Cluster compactness and separation

**Visualization**:
- **2D Projections**: PCA, t-SNE scatter plots
- **Cluster Maps**: Color-coded cluster assignments
- **Dendrograms**: Hierarchical clustering trees
- **Feature Importance**: Contribution of each descriptor

### 14.6 Model Selection

**Cluster Number Optimization**:
1. **Elbow Method**: Plot inertia vs. number of clusters
2. **Silhouette Analysis**: Evaluate cluster quality scores
3. **Gap Statistic**: Compare with random data clustering
4. **Chemical Validation**: Ensure chemically meaningful groups

**Cross-Validation**:
- **Stability**: Consistency across data subsets
- **Robustness**: Sensitivity to parameter changes
- **Reproducibility**: Multiple random initializations

### 14.7 Filtering and Selection

**Cluster-Based Filtering**:
- **Representative Selection**: Choose cluster centroids
- **Diversity Sampling**: Select from each cluster
- **Property-Based**: Filter by thermodynamic properties
- **Size-Based**: Consider cluster populations

**Interactive Selection**:
- **Manual Curation**: Expert-guided selection
- **Property Constraints**: Apply chemical filters
- **Visualization-Guided**: Use plot interactions
- **Iterative Refinement**: Multiple selection rounds

### 14.8 Workflow

1. **Load Feature Data**: Import molecular descriptors from previous steps
2. **Preprocess Features**: Normalize and select relevant descriptors
3. **Choose Algorithm**: Select clustering method and parameters
4. **Perform Clustering**: Execute clustering analysis
5. **Evaluate Results**: Assess cluster quality and chemical meaning
6. **Visualize Clusters**: Generate plots and projections
7. **Select Representatives**: Choose molecules for further study
8. **Export Results**: Save cluster assignments and selections

### 14.9 Output Files

**Cluster Assignments** (`clusters_*.csv`):
- Molecule-to-cluster mapping
- Cluster centroids and properties
- Quality metrics and statistics

**Visualization Files** (`clusters_*.svg`):
- 2D cluster projections
- Dendrogram plots
- Feature importance charts

**Selection Results** (`sel_model_indices.pomsim`):
- Selected molecule indices
- Selection criteria and rationale
- Representative structures per cluster

**Filtered Data** (`filtering/filt_by_*.svg`):
- Cluster-specific analysis plots
- Property distributions per cluster
- Comparative visualizations

---

## 15. Plotting Tab

### 15.1 Overview

> **Figure 11.** Plotting Tab Interface

The Plotting tab creates publication-quality visualizations of simulation results. It provides comprehensive plotting capabilities for speciation diagrams, phase diagrams, and statistical analyses.

### 15.2 Plot Types

**Speciation Plots**:
- **Alpha Diagrams**: Species fraction vs. pH
- **Concentration Plots**: Absolute concentrations vs. conditions
- **3D Surfaces**: Multi-variable speciation landscapes
- **Contour Maps**: 2D projections of 3D data

**Phase Diagrams**:
- **Predominance Diagrams**: Dominant species regions
- **Stability Fields**: Thermodynamic stability boundaries
- **Pourbaix Diagrams**: pH-potential relationships
- **Temperature-Composition**: Phase behavior vs. temperature

**Statistical Plots**:
- **Correlation Plots**: Computed vs. experimental values
- **Residual Analysis**: Error distribution patterns
- **Cluster Visualizations**: Machine learning results
- **Feature Importance**: Descriptor contribution analysis

### 15.3 Customization Options

**Appearance Settings**:

**Colors and Styles**:
- **Color Schemes**: Predefined scientific color palettes
- **Line Styles**: Solid, dashed, dotted patterns
- **Markers**: Various point styles and sizes
- **Transparency**: Alpha channel control

**Axes and Labels**:
- **Axis Ranges**: Manual or automatic scaling
- **Tick Marks**: Spacing and formatting
- **Labels**: Font size, style, and positioning
- **Units**: Automatic unit conversion and display

**Layout Options**:
- **Figure Size**: Width and height in inches or pixels
- **Margins**: Border spacing around plots
- **Subplots**: Multiple panels in single figure
- **Legends**: Position, style, and content

### 15.4 Export Formats

**Raster Formats**:
- **PNG**: High-quality bitmap images
- **JPEG**: Compressed photographs
- **TIFF**: Uncompressed scientific images
- **Resolution**: Customizable DPI settings

**Vector Formats**:
- **SVG**: Scalable vector graphics
- **PDF**: Publication-ready documents
- **EPS**: Encapsulated PostScript
- **EMF**: Enhanced metafile (Windows)

**Data Formats**:
- **CSV**: Numerical data tables
- **JSON**: Structured data export
- **HDF5**: Large dataset storage
- **Excel**: Spreadsheet-compatible format

### 15.5 Interactive Features

**Zoom and Pan**:
- **Mouse Controls**: Wheel zoom, drag pan
- **Keyboard Shortcuts**: Arrow keys for fine control
- **Reset View**: Return to original scale
- **Fit to Data**: Automatic optimal scaling

**Data Inspection**:
- **Hover Information**: Point values on mouse over
- **Click Selection**: Highlight specific data points
- **Crosshairs**: Precise coordinate reading
- **Measurement Tools**: Distance and angle measurement

### 15.6 Batch Processing

**Multiple Plots**:
- **Series Generation**: Automatic plot series creation
- **Parameter Sweeps**: Vary conditions systematically
- **Comparison Plots**: Side-by-side visualizations
- **Animation**: Time-series or parameter evolution

**Template System**:
- **Style Templates**: Predefined formatting schemes
- **Layout Templates**: Standard figure arrangements
- **Custom Templates**: User-defined plot styles
- **Template Sharing**: Export/import formatting

### 15.7 Quality Control

**Publication Standards**:
- **Resolution**: Minimum 300 DPI for print
- **Font Sizes**: Readable at publication scale
- **Color Blindness**: Accessible color schemes
- **Contrast**: Sufficient for black/white printing

**Validation Checks**:
- **Data Integrity**: Verify plot data accuracy
- **Label Completeness**: Ensure all elements labeled
- **Unit Consistency**: Check dimensional analysis
- **Scale Appropriateness**: Verify axis ranges

### 15.8 Workflow

1. **Select Data Source**: Choose results from previous analysis steps
2. **Choose Plot Type**: Select appropriate visualization method
3. **Configure Appearance**: Set colors, styles, and layout
4. **Preview Plot**: Review appearance and make adjustments
5. **Export Plot**: Save in desired format and resolution
6. **Batch Processing**: Generate multiple related plots
7. **Quality Check**: Verify publication readiness

### 15.9 Output Files

**Plot Images** (`*.png`, `*.svg`, `*.pdf`):
- High-resolution publication-quality figures
- Multiple format options for different uses
- Consistent styling and professional appearance

**Plot Data** (`plot_data_*.csv`):
- Numerical data underlying each plot
- Reproducibility and further analysis
- Custom plot generation capability

**Style Sheets** (`plot_styles_*.json`):
- Formatting parameters and preferences
- Template sharing and consistency
- Automated styling application

---

## 16. Configuration Management

### 16.1 Overview

Configuration management allows you to save and restore complete simulation setups, ensuring reproducibility and enabling parameter sharing between projects.

### 16.2 Configuration Files

**File Format**: `.pomsim` files use standard INI format with sections for each tab
**Content**: All GUI parameters, file paths, and calculation settings
**Portability**: Configurations can be shared between users and systems
**Versioning**: Files include version information for compatibility checking

### 16.3 Saving Configurations

**Export Process**:
1. **Complete Setup**: Ensure all desired parameters are configured
2. **File Menu**: Select "Export Configuration File..."
3. **Choose Location**: Navigate to desired save directory
4. **Name File**: Use descriptive filename with `.pomsim` extension
5. **Confirm Save**: Verify successful file creation

**What Gets Saved**:
- All input field values across all tabs
- File path selections and directory settings
- Algorithm choices and parameter values
- Display preferences and formatting options

### 16.4 Loading Configurations

**Import Process**:
1. **File Menu**: Select "Import Configuration File..."
2. **Browse Files**: Navigate to configuration file location
3. **Select File**: Choose `.pomsim` file to load
4. **Automatic Population**: GUI fields are automatically filled
5. **Verify Settings**: Review loaded parameters for accuracy

**Validation Checks**:
- **File Paths**: Verify referenced files exist
- **Parameter Ranges**: Check values are within valid limits
- **Compatibility**: Ensure version compatibility
- **Dependencies**: Validate required data availability

### 16.5 Sample Configurations

The application includes several sample configurations in the `inputs` folder:

**config_W.pomsim**: Tungsten polyoxometalate system
- Complete workflow setup for tungsten species
- Experimental data from literature sources
- Optimized parameters for W-based systems

**config_PMo.pomsim**: Phosphomolybdate system
- Heteropolyanion configuration example
- Multi-metal system parameters
- Advanced clustering and analysis settings

**config_As.pomsim**: Arsenate system
- Simple isopolyanion example
- Basic workflow demonstration
- Educational and testing purposes

**config_C.pomsim**: Carbonate system
- Alternative chemistry example
- Different calculation approaches
- Comparative analysis setup

### 16.6 Best Practices

**File Organization**:
- Use descriptive filenames indicating system and purpose
- Organize configurations by project or chemical system
- Include date and version information in filenames
- Maintain backup copies of important configurations

**Documentation**:
- Add comments to configuration files when possible
- Keep notes about parameter choices and rationale
- Document data sources and experimental conditions
- Record any modifications from standard procedures

**Sharing and Collaboration**:
- Include all necessary data files with configurations
- Use relative paths when possible for portability
- Provide clear instructions for configuration use
- Validate configurations on different systems

### 16.7 Troubleshooting

**Common Issues**:

**Missing Files**:
- **Problem**: Referenced files not found
- **Solution**: Update file paths or copy missing files
- **Prevention**: Use relative paths and include all dependencies

**Parameter Conflicts**:
- **Problem**: Incompatible parameter combinations
- **Solution**: Review and adjust conflicting settings
- **Prevention**: Validate configurations before saving

**Version Incompatibility**:
- **Problem**: Configuration from different software version
- **Solution**: Update parameters to current format
- **Prevention**: Include version information in filenames

---

## 17. Console Output

### 17.1 Overview

> **Figure 12.** Console Output Panel

The console panel provides real-time feedback about application operations, displaying progress updates, warnings, and error messages.

### 17.2 Message Types

**Information Messages**:
- **Format**: Standard black text
- **Content**: Progress updates, successful operations
- **Examples**: "Loading configuration file...", "Calculation completed successfully"

**Warning Messages**:
- **Format**: Orange or yellow text
- **Content**: Non-critical issues, recommendations
- **Examples**: "Parameter outside recommended range", "Missing optional data"

**Error Messages**:
- **Format**: Red text
- **Content**: Critical problems, failed operations
- **Examples**: "File not found", "Calculation failed due to convergence issues"

**Debug Messages**:
- **Format**: Gray text
- **Content**: Detailed technical information
- **Visibility**: Can be enabled/disabled in preferences

### 17.3 Console Features

**Scrolling**: Automatic scrolling to show latest messages
**Search**: Find specific messages using Ctrl+F
**Copy**: Select and copy messages for reporting
**Clear**: Clear console history with dedicated button
**Export**: Save console log to text file

### 17.4 Progress Monitoring

**Operation Status**:
- Real-time updates during long calculations
- Percentage completion where applicable
- Estimated time remaining for lengthy operations
- Step-by-step progress through complex workflows

**Resource Usage**:
- Memory consumption warnings
- CPU utilization information
- Disk space requirements
- Network activity for data downloads

### 17.5 Error Diagnosis

**Error Categories**:

**File System Errors**:
- Missing input files
- Permission issues
- Disk space problems
- Path resolution failures

**Calculation Errors**:
- Convergence failures
- Numerical instabilities
- Parameter out of range
- Insufficient data

**System Errors**:
- Memory allocation failures
- Library import problems
- Version compatibility issues
- Hardware limitations

### 17.6 Console Management

**Size Control**: Resize console panel by dragging splitter
**Theme Integration**: Console colors match application theme
**Font Settings**: Monospace font for better readability
**History Limit**: Configurable maximum number of stored messages

---

## 18. Keyboard Shortcuts

### 18.1 File Operations

| Shortcut | Action | Description |
|----------|--------|-------------|
| `Ctrl+O` | Import Configuration | Open configuration file dialog |
| `Ctrl+S` | Export Configuration | Save current settings to file |
| `Ctrl+Q` | Exit Application | Close POMSimulator |

### 18.2 Interface Control

| Shortcut | Action | Description |
|----------|--------|-------------|
| `Ctrl+T` | Toggle Theme | Switch between light and dark themes |
| `Ctrl+M` | Toggle Mode | Switch between simulation and visualization modes |
| `Ctrl+F5` | Refresh GUI | Reset interface to default state |
| `Ctrl+P` | Preferences | Open application settings |

### 18.3 Help and Documentation

| Shortcut | Action | Description |
|----------|--------|-------------|
| `F1` | Context Help | Show help for current tab |
| `F2` | Documentation | Open online documentation |
| `Ctrl+F1` | Key Bindings | Display keyboard shortcuts |
| `Ctrl+A` | About | Show application information |

### 18.4 Text Editing

| Shortcut | Action | Description |
|----------|--------|-------------|
| `Ctrl+C` | Copy | Copy selected text |
| `Ctrl+V` | Paste | Paste from clipboard |
| `Ctrl+X` | Cut | Cut selected text |
| `Ctrl+Z` | Undo | Undo last action |
| `Ctrl+Y` | Redo | Redo last undone action |

### 18.5 Navigation

| Shortcut | Action | Description |
|----------|--------|-------------|
| `Tab` | Next Field | Move to next input field |
| `Shift+Tab` | Previous Field | Move to previous input field |
| `Enter` | Activate | Activate focused button or field |
| `Escape` | Cancel | Cancel current dialog or operation |

### 18.6 Console Operations

| Shortcut | Action | Description |
|----------|--------|-------------|
| `Ctrl+F` | Find in Console | Search console messages |
| `Ctrl+L` | Clear Console | Clear all console messages |
| `Ctrl+Shift+C` | Copy Console | Copy all console text |

---

## 19. Troubleshooting

### 19.1 Installation Issues

**Python Environment Problems**:
- **Symptom**: Import errors or missing modules
- **Solution**: Verify Python version (3.7+) and install required packages
- **Command**: `pip install -r requirements.txt`

**PyQt5 Installation**:
- **Symptom**: GUI fails to start or displays incorrectly
- **Solution**: Reinstall PyQt5 with proper system compatibility
- **Command**: `pip uninstall PyQt5 && pip install PyQt5`

**Permission Issues**:
- **Symptom**: Cannot write output files or access directories
- **Solution**: Run with appropriate permissions or change output directory
- **Prevention**: Use user-writable directories for outputs

### 19.2 File System Issues

**Missing Input Files**:
- **Symptom**: "File not found" errors during operations
- **Solution**: Verify file paths and ensure all required files are present
- **Check**: Use file browser to confirm file existence and accessibility

**Path Resolution Problems**:
- **Symptom**: Incorrect file paths or directory access failures
- **Solution**: Use absolute paths or verify relative path correctness
- **Tip**: Avoid special characters and spaces in file paths

**Large File Handling**:
- **Symptom**: Memory errors or slow performance with large datasets
- **Solution**: Increase available memory or process data in smaller chunks
- **Monitoring**: Watch memory usage in system task manager

### 19.3 Calculation Errors

**Convergence Failures**:
- **Symptom**: Calculations fail to reach solution
- **Solution**: Adjust convergence criteria or initial guesses
- **Parameters**: Increase iteration limits or change tolerance values

**Numerical Instabilities**:
- **Symptom**: Unrealistic results or calculation crashes
- **Solution**: Check input data quality and parameter ranges
- **Validation**: Compare results with known experimental values

**Memory Limitations**:
- **Symptom**: Out of memory errors during large calculations
- **Solution**: Reduce problem size or increase available RAM
- **Optimization**: Close other applications to free memory

### 19.4 Display Issues

**Scaling Problems**:
- **Symptom**: Interface elements too small or large
- **Solution**: Adjust system DPI settings or use application scaling options
- **Settings**: Check display scaling in system preferences

**Theme Issues**:
- **Symptom**: Poor contrast or unreadable text
- **Solution**: Switch themes or adjust system color settings
- **Shortcut**: Use `Ctrl+T` to toggle between light and dark themes

**Plot Display Problems**:
- **Symptom**: Plots not displaying or appearing corrupted
- **Solution**: Update graphics drivers or change plot backend
- **Alternative**: Export plots to files for external viewing

### 19.5 Performance Issues

**Slow Startup**:
- **Symptom**: Application takes long time to launch
- **Solution**: Check for antivirus interference or disk space issues
- **Optimization**: Close unnecessary background applications

**Slow Calculations**:
- **Symptom**: Operations take much longer than expected
- **Solution**: Verify system meets minimum requirements
- **Monitoring**: Check CPU and memory usage during operations

**Unresponsive Interface**:
- **Symptom**: GUI freezes during operations
- **Solution**: Use Stop buttons to cancel long-running operations
- **Prevention**: Monitor progress and avoid extremely large calculations

### 19.6 Data Issues

**Corrupted Files**:
- **Symptom**: Unexpected errors when loading data files
- **Solution**: Verify file integrity and re-generate if necessary
- **Backup**: Maintain backup copies of important data

**Format Incompatibilities**:
- **Symptom**: Files cannot be read or produce incorrect results
- **Solution**: Check file format specifications and conversion tools
- **Validation**: Use sample files to verify format compatibility

**Missing Dependencies**:
- **Symptom**: Features unavailable or import errors
- **Solution**: Install optional dependencies for full functionality
- **Example**: ASE for molecular visualization, additional plotting libraries

---

## 20. Frequently Asked Questions

### 20.1 General Usage

**Q: What types of chemical systems can POMSimulator handle?**
A: POMSimulator is designed for polyoxometalate systems, including both isopolyanions (single metal type) and heteropolyanions (multiple metal types). It supports tungsten, molybdenum, phosphorus, arsenic, and carbon-based systems.

**Q: Do I need quantum chemistry experience to use POMSimulator?**
A: Basic understanding of quantum chemistry concepts is helpful, but the software provides guided workflows and sample data for learning. The interface is designed to be accessible to users with varying levels of expertise.

**Q: Can I use experimental data from different literature sources?**
A: Yes, the software allows you to input custom experimental formation constants. Ensure data consistency in terms of ionic strength, temperature, and standard states for best results.

### 20.2 Technical Questions

**Q: What file formats does POMSimulator support?**
A: The software supports ADF output files (.out), MOL molecular structure files (.mol), configuration files (.pomsim), and various data formats (CSV, TXT, NPZ). It can export plots in PNG, SVG, PDF, and other standard formats.

**Q: How accurate are the calculated formation constants?**
A: Accuracy depends on the quality of input quantum chemistry data and the scaling procedure. Properly calibrated results typically achieve agreement within 1-2 log units of experimental values for well-studied systems.

**Q: Can I run calculations on multiple systems simultaneously?**
A: The current version processes one system at a time within the GUI. However, you can run multiple instances of the application or use batch processing scripts for parallel calculations.

### 20.3 Workflow Questions

**Q: Do I need to complete all tabs in order?**
A: The workflow is designed to be sequential (Presimulation → Simulation → Scaling → Speciation/Clustering → Plotting), but you can skip certain steps if you have pre-existing data. The Clustering and Speciation tabs can be run independently after Scaling.

**Q: How do I know if my calculations are successful?**
A: Check the console output for completion messages, verify that output files are generated in the expected locations, and review results for chemical reasonableness. The interface provides visual indicators of completion status.

**Q: Can I modify parameters and re-run calculations?**
A: Yes, you can adjust parameters at any stage and re-run calculations. The software will overwrite previous results, so save important configurations and results before making changes.

### 20.4 Data Management

**Q: Where are my results saved?**
A: Results are saved in the `outputs` folder within the project directory, organized by system name (e.g., `W_data`, `PMo_data`). Each system has its own subfolder containing all generated files.

**Q: How do I share my work with collaborators?**
A: Export your configuration using the File menu to create a `.pomsim` file containing all parameters. Share this file along with any custom input data files. Collaborators can import the configuration to reproduce your setup.

**Q: Can I export data for use in other software?**
A: Yes, most numerical results are saved in CSV format for easy import into spreadsheet applications, statistical software, or custom analysis scripts. Plots can be exported in various formats for publication use.

### 20.5 Troubleshooting

**Q: The application won't start. What should I do?**
A: First, verify that Python 3.7+ is installed and all required packages are available. Check the console for error messages. Try running from the command line to see detailed error information.

**Q: My calculations are taking a very long time. Is this normal?**
A: Calculation time depends on system size and complexity. Large systems with many species can take several minutes to hours. Monitor the console for progress updates and use the Stop button if needed.

**Q: The interface looks wrong or is difficult to read.**
A: Try toggling the theme using `Ctrl+T` or adjusting your system's display scaling settings. The application automatically adapts to high-DPI displays, but manual adjustment may be needed in some cases.

### 20.6 Advanced Usage

**Q: Can I add my own experimental data to the database?**
A: Yes, you can modify the experimental constants in the database files or input custom values through the interface. Ensure proper formatting and units for consistency.

**Q: How do I cite POMSimulator in publications?**
A: Citation information is available in the About dialog (`Ctrl+A`). Include the software version, authors, and any relevant publications describing the methods implemented.

**Q: Are there plans for additional features?**
A: The software is actively developed with new features added based on user feedback and scientific needs. Check the documentation and project repository for updates and roadmap information.

---

## Appendices

### Appendix A: Sample Data Description

The `inputs` folder contains several example datasets:

- **W_Set_PBE**: Tungsten polyoxometalate quantum chemistry outputs
- **PMo_Set**: Phosphomolybdate heteropolyanion data
- **As_Set**: Arsenate isopolyanion examples
- **C_Set**: Carbonate system for comparison

### Appendix B: File Format Specifications

Detailed specifications for supported file formats are available in the technical documentation.

### Appendix C: Algorithm References

Scientific references for the computational methods implemented in POMSimulator.

### Appendix D: Validation Studies

Results from validation studies comparing POMSimulator predictions with experimental data.

---

**End of Manual**

*This manual provides comprehensive guidance for using POMSimulator. For additional support, consult the online documentation or contact the development team.*