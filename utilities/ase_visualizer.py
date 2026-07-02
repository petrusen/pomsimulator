"""
Molecule Visualizer using ASE

This module provides functionality to visualize .mol files using the Atomic Simulation Environment (ASE).
"""

import os
import sys
import numpy as np
from ase import Atoms
from ase.io import read
from ase.visualize import view
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
                            QComboBox, QFileDialog, QCheckBox, QSlider, QColorDialog, QMessageBox)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFontMetrics
from pathlib import Path


class ElidedLabel(QLabel):
    """A QLabel that automatically elides text when it doesn't fit."""
    
    def __init__(self, text="", parent=None):
        super().__init__(text, parent)
        self._full_text = text
        
    def setText(self, text):
        """Set the text and store the full version for elision."""
        self._full_text = text
        self._update_elided_text()
        
    def _update_elided_text(self):
        """Update the displayed text with elision if necessary."""
        if not self._full_text:
            super().setText("")
            return
            
        # Get available width
        available_width = self.width() - 10  # Leave some margin
        if available_width <= 0:
            available_width = 200  # Default reasonable width
            
        # Use font metrics to elide text
        font_metrics = QFontMetrics(self.font())
        elided_text = font_metrics.elidedText(self._full_text, Qt.ElideMiddle, available_width)
        super().setText(elided_text)
        
    def resizeEvent(self, event):
        """Handle resize events to update elided text."""
        super().resizeEvent(event)
        self._update_elided_text()
        
    def get_full_text(self):
        """Get the full, non-elided text."""
        return self._full_text

class MoleculeCanvas(FigureCanvas):
    """
    A canvas for displaying molecules using ASE and Matplotlib.
    """
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111, projection='3d')
        
        super(MoleculeCanvas, self).__init__(self.fig)
        self.setParent(parent)
        
        self.molecule = None
        self.atom_colors = {
            'H': 'white',
            'C': 'black',
            'N': 'blue',
            'O': 'red',
            'P': 'orange',
            'S': 'yellow',
            'Mo': 'green',
            'W': 'purple',
            'V': 'cyan'
        }
        
    def load_molecule(self, mol_file):
        """Load a molecule from a .mol file"""
        try:
            # Read the molecule using ASE
            self.molecule = read(mol_file, format='mol')
            self.plot_molecule()
            return True
        except Exception as e:
            print(f"Error loading molecule: {str(e)}")
            return False
            
    def plot_molecule(self):
        """Plot the molecule in 3D"""
        if self.molecule is None:
            return
            
        self.axes.clear()
        
        # Get positions and elements
        positions = self.molecule.get_positions()
        elements = self.molecule.get_chemical_symbols()
        
        # Plot atoms
        for i, (pos, element) in enumerate(zip(positions, elements)):
            color = self.atom_colors.get(element, 'gray')
            size = 200 if element in ['Mo', 'W', 'V'] else 100
            self.axes.scatter(pos[0], pos[1], pos[2], c=color, s=size, edgecolors='black')
            
        # Plot bonds
        bonds = self.get_bonds_from_molecule()
        for i, j in bonds:
            pos_i = positions[i]
            pos_j = positions[j]
            self.axes.plot([pos_i[0], pos_j[0]], 
                          [pos_i[1], pos_j[1]], 
                          [pos_i[2], pos_j[2]], 'k-', lw=1)
        
        # Set equal aspect ratio
        self.axes.set_box_aspect([1, 1, 1])
        
        # Remove axis labels and ticks
        self.axes.set_axis_off()
        
        self.fig.tight_layout()
        self.draw()
        
    def get_bonds_from_molecule(self):
        """Extract bonds from the molecule"""
        bonds = []
        
        # If the molecule has connectivity information
        if hasattr(self.molecule, 'get_all_distances'):
            # Get all pairwise distances
            distances = self.molecule.get_all_distances()
            elements = self.molecule.get_chemical_symbols()
            
            # Define bond length thresholds based on elements
            bond_thresholds = {
                ('H', 'H'): 0.8,
                ('C', 'C'): 1.6,
                ('C', 'H'): 1.2,
                ('C', 'O'): 1.5,
                ('O', 'H'): 1.1,
                ('O', 'O'): 1.5,
                ('Mo', 'O'): 2.2,
                ('W', 'O'): 2.2,
                ('V', 'O'): 2.0,
            }
            
            # Default threshold for other element pairs
            default_threshold = 2.0
            
            # Find bonds based on distance thresholds
            n_atoms = len(self.molecule)
            for i in range(n_atoms):
                for j in range(i+1, n_atoms):
                    el_i = elements[i]
                    el_j = elements[j]
                    
                    # Get threshold for this element pair
                    pair = tuple(sorted([el_i, el_j]))
                    threshold = bond_thresholds.get(pair, default_threshold)
                    
                    # If distance is less than threshold, add bond
                    if distances[i, j] < threshold:
                        bonds.append((i, j))
        
        return bonds
        
    def rotate_view(self, azimuth, elevation):
        """Rotate the view of the molecule"""
        if self.molecule is None:
            return
            
        self.axes.view_init(elev=elevation, azim=azimuth)
        self.draw()
        
    def set_atom_color(self, element, color):
        """Set the color for a specific element"""
        self.atom_colors[element] = color
        self.plot_molecule()


class MoleculeVisualizer(QWidget):
    """
    A widget for visualizing molecules from .mol files.
    """
    def __init__(self, parent=None):
        super(MoleculeVisualizer, self).__init__(parent)
        self._full_file_path = None  # Store the complete file path internally
        self.init_ui()
        
    def init_ui(self):
        """
        Initialize the user interface for the molecule visualizer.
        
        This method sets up the complete UI layout for the molecule visualizer, including:
        - The 3D molecule canvas for visualization
        - File selection controls for loading molecule files
        - View controls with sliders for adjusting azimuth and elevation angles
        - Element color controls for customizing the appearance of different atoms
        
        The UI is organized in a vertical layout with the molecule canvas at the top
        and control panels at the bottom. The control panels are arranged horizontally
        and include file selection, view controls, and element color customization.
        
        Parameters:
            None
            
        Returns:
            None: This method sets up the UI components but doesn't return a value.
        """
        layout = QVBoxLayout()
        
        # Molecule canvas
        self.canvas = MoleculeCanvas(self, width=8, height=6)
        layout.addWidget(self.canvas)
        
        # Controls
        controls_layout = QHBoxLayout()
        
        # File selection
        file_layout = QVBoxLayout()
        file_label = QLabel("Molecule File:")
        file_layout.addWidget(file_label)
        
        file_select_layout = QHBoxLayout()
        self.file_path = ElidedLabel("No file selected")
        # Configure the label to handle long paths gracefully
        self.file_path.setWordWrap(False)
        self.file_path.setTextInteractionFlags(Qt.TextSelectableByMouse)
        # Set a reasonable minimum width and allow shrinking
        self.file_path.setMinimumWidth(100)
        self.file_path.setMaximumWidth(300)  # Prevent excessive expansion
        self.file_path.setSizePolicy(self.file_path.sizePolicy().horizontalPolicy(),
                                   self.file_path.sizePolicy().verticalPolicy())
        file_select_layout.addWidget(self.file_path)
        
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(self.browse_file)
        file_select_layout.addWidget(browse_btn)
        
        file_layout.addLayout(file_select_layout)
        controls_layout.addLayout(file_layout)
        
        # View controls
        view_layout = QVBoxLayout()
        view_label = QLabel("View Controls:")
        view_layout.addWidget(view_label)
        
        # Azimuth slider
        azimuth_layout = QHBoxLayout()
        azimuth_layout.addWidget(QLabel("Azimuth:"))
        self.azimuth_slider = QSlider(Qt.Horizontal)
        self.azimuth_slider.setRange(0, 360)
        self.azimuth_slider.setValue(30)
        self.azimuth_slider.valueChanged.connect(self.update_view)
        azimuth_layout.addWidget(self.azimuth_slider)
        view_layout.addLayout(azimuth_layout)
        
        # Elevation slider
        elevation_layout = QHBoxLayout()
        elevation_layout.addWidget(QLabel("Elevation:"))
        self.elevation_slider = QSlider(Qt.Horizontal)
        self.elevation_slider.setRange(-90, 90)
        self.elevation_slider.setValue(30)
        self.elevation_slider.valueChanged.connect(self.update_view)
        elevation_layout.addWidget(self.elevation_slider)
        view_layout.addLayout(elevation_layout)
        
        controls_layout.addLayout(view_layout)
        
        # Element color controls
        color_layout = QVBoxLayout()
        color_label = QLabel("Element Colors:")
        color_layout.addWidget(color_label)
        
        element_layout = QHBoxLayout()
        element_layout.addWidget(QLabel("Element:"))
        self.element_combo = QComboBox()
        self.element_combo.addItems(['H', 'C', 'N', 'O', 'P', 'S', 'Mo', 'W', 'V'])
        element_layout.addWidget(self.element_combo)
        
        color_btn = QPushButton("Change Color")
        color_btn.clicked.connect(self.change_element_color)
        element_layout.addWidget(color_btn)
        
        color_layout.addLayout(element_layout)
        controls_layout.addLayout(color_layout)
        
        # Add close button
        close_layout = QVBoxLayout()
        close_label = QLabel("Actions:")
        close_layout.addWidget(close_label)
        
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close_viewer)
        close_layout.addWidget(close_btn)
        
        controls_layout.addLayout(close_layout)
        
        layout.addLayout(controls_layout)
        
        self.setLayout(layout)
        self.setWindowTitle("Molecule Visualizer")
        
    def browse_file(self):
        """Open a file dialog to select a .mol file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Molecule File", "", "Molecule Files (*.mol)",
        options=QFileDialog.DontUseNativeDialog)
        
        if file_path:
            self.set_file_path(file_path)
            success = self.canvas.load_molecule(file_path)
            if success:
                self.update_view()
            else:
                self.file_path.setText("Error loading file")
                self.file_path.setToolTip("Failed to load the selected molecule file")
                self._full_file_path = None
                
    def update_view(self):
        """Update the view based on slider values"""
        azimuth = self.azimuth_slider.value()
        elevation = self.elevation_slider.value()
        self.canvas.rotate_view(azimuth, elevation)
        
    def change_element_color(self):
        """Open a color dialog to change element color"""
        element = self.element_combo.currentText()
        current_color = self.canvas.atom_colors.get(element, 'gray')
        
        color = QColorDialog.getColor()
        if color.isValid():
            self.canvas.set_atom_color(element, color.name())
            
    def close_viewer(self):
        """Close the molecule viewer."""
        # Find the PreviewManager parent and close the preview
        parent = self.parent()
        while parent is not None:
            # Import here to avoid circular imports
            from pomsimulator.GUI import PreviewManager
            if isinstance(parent, PreviewManager):
                parent.close_preview()
                return
            parent = parent.parent()
        
        # Fallback: hide this widget directly
        self.hide()
        
    def set_file_path(self, file_path):
        """
        Set the file path, storing the full path internally and displaying a shortened version.
        
        Parameters:
            file_path (str): The full path to the molecule file.
        """
        if not file_path:
            self._full_file_path = None
            self.file_path.setText("No file selected")
            self.file_path.setToolTip("")
            return
            
        # Store the full path internally
        self._full_file_path = file_path
        
        # Create a shortened display version
        path_obj = Path(file_path)
        
        # Try to show parent directory + filename for context, or just filename if too long
        try:
            parent_name = path_obj.parent.name
            filename = path_obj.name
            
            if parent_name and len(parent_name) < 20:  # Reasonable parent directory name length
                display_text = f"{parent_name}/{filename}"
            else:
                display_text = filename
        except:
            # Fallback to just filename if path operations fail
            display_text = path_obj.name if hasattr(path_obj, 'name') else os.path.basename(file_path)
        
        # Set the display text and tooltip
        self.file_path.setText(display_text)
        self.file_path.setToolTip(f"Full path: {file_path}")
        
    def get_file_path(self):
        """
        Get the full file path.
        
        Returns:
            str: The complete file path, or None if no file is loaded.
        """
        return self._full_file_path

    def open_mol_file(self, file_path):
        """
        Open a .mol file in the molecule visualizer.

        This method launches the molecule visualizer and loads the specified .mol file.

        Parameters:
            file_path (str): The full path to the .mol file to be opened.

        Returns:
            None: This method displays a visualizer window but doesn't return a value.
        """
        try:
            from pomsimulator.utilities.ase_visualizer import MoleculeVisualizer

            # Create and show the visualizer
            self.visualizer = MoleculeVisualizer()
            
            # Load the molecule file
            self.visualizer.file_path.setText(file_path)
            success = self.visualizer.canvas.load_molecule(file_path)

            if not success:
                QMessageBox.warning(self, "Error Loading Molecule",
                                   "Could not load the molecule file. It may be in an unsupported format.")
            else:
                # Only show the visualizer if the molecule was loaded successfully
                self.visualizer.show()

        except ImportError as e:
            QMessageBox.warning(self, "Import Error",
                               f"Could not load the molecule visualizer: {str(e)}\n\n"
                               "Make sure ASE (Atomic Simulation Environment) is installed.")


def parse_mol_file(mol_file):
    """
    Parse a .mol file and return atom positions and bonds.
    
    Args:
        mol_file (str): Path to the .mol file
        
    Returns:
        tuple: (atoms, bonds) where atoms is a list of (element, x, y, z) tuples
               and bonds is a list of (atom1_idx, atom2_idx, bond_type) tuples
    """
    with open(mol_file, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    counts_line = lines[3]
    n_atoms = int(counts_line[0:3].strip())
    n_bonds = int(counts_line[3:6].strip())
    
    # Parse atoms
    atoms = []
    for i in range(n_atoms):
        line = lines[4 + i]
        x = float(line[0:10].strip())
        y = float(line[10:20].strip())
        z = float(line[20:30].strip())
        element = line[31:34].strip()
        atoms.append((element, x, y, z))
    
    # Parse bonds
    bonds = []
    for i in range(n_bonds):
        line = lines[4 + n_atoms + i]
        atom1 = int(line[0:3].strip()) - 1  # Mol files are 1-indexed
        atom2 = int(line[3:6].strip()) - 1
        bond_type = int(line[6:9].strip())
        bonds.append((atom1, atom2, bond_type))
    
    return atoms, bonds


def create_ase_atoms_from_mol(mol_file):
    """
    Create an ASE Atoms object from a .mol file.
    
    Args:
        mol_file (str): Path to the .mol file
        
    Returns:
        ase.Atoms: ASE Atoms object representing the molecule
    """
    atoms, bonds = parse_mol_file(mol_file)
    
    # Extract elements and positions
    elements = [atom[0] for atom in atoms]
    positions = np.array([[atom[1], atom[2], atom[3]] for atom in atoms])
    
    # Create ASE Atoms object
    atoms_obj = Atoms(symbols=elements, positions=positions)
    
    return atoms_obj


def launch_visualizer(mol_file=None, standalone=True):
    """
    Launch the molecule visualizer.
    
    Args:
        mol_file (str, optional): Path to a .mol file to load initially
        standalone (bool): Whether to run as a standalone app (True) or return the visualizer (False)
    
    Returns:
        If standalone is True, doesn't return (calls sys.exit())
        If standalone is False, returns the MoleculeVisualizer instance
    """
    from PyQt5.QtWidgets import QApplication
    
    # Only create a QApplication if running standalone or if one doesn't exist
    if standalone and not QApplication.instance():
        # Enable High-DPI support before creating QApplication
        from pomsimulator.dpi_utils import enable_high_dpi_support
        enable_high_dpi_support()
        app = QApplication(sys.argv)
    
    visualizer = MoleculeVisualizer()
    
    # Load molecule if provided
    if mol_file and os.path.isfile(mol_file):
        visualizer.file_path.setText(mol_file)
        visualizer.canvas.load_molecule(mol_file)
        visualizer.update_view()
    
    if standalone:
        visualizer.show()
        sys.exit(app.exec_())
    else:
        return visualizer


if __name__ == "__main__":
    # Only run this code if the script is executed directly
    from PyQt5.QtWidgets import QApplication
    import sys
    
    # Enable High-DPI support before creating QApplication
    from pomsimulator.dpi_utils import enable_high_dpi_support
    enable_high_dpi_support()
    
    app = QApplication(sys.argv)
    visualizer = MoleculeVisualizer()
    visualizer.show()
    sys.exit(app.exec_())
