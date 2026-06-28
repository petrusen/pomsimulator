import os
import subprocess
import sys
from configparser import ConfigParser
from multiprocessing import cpu_count
from io import StringIO
import traceback
import tempfile
from datetime import datetime

from PyQt5.QtCore import Qt, QModelIndex, QLocale, QThread, pyqtSignal, QTimer, QUrl
from PyQt5.QtGui import QFont, QIcon, QPixmap, QDesktopServices, QValidator, QPalette, QColor, QKeySequence
from PyQt5.QtWidgets import (QApplication, QMainWindow, QTabWidget, QWidget, QVBoxLayout,
                             QHBoxLayout, QGridLayout, QLineEdit, QComboBox, QFileDialog, QCheckBox, QSpinBox,
                             QDoubleSpinBox, QGroupBox, QRadioButton, QMessageBox, QTreeView, QFileSystemModel,
                             QDockWidget, QDialog, QToolBar, QButtonGroup, QPlainTextEdit, QSlider, QProgressBar,
                             QTextEdit, QLabel, QAction, QPushButton, QScrollArea, QShortcut, QFrame)

from pomsimulator.modules.DataBase import experimental_constants, allowed_scaling_modes, reaction_references, \
    clustering_features

# Set the script directory for relative imports
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_dir)


class POMSimulatorGUI(QMainWindow):

    def __init__(self):
        """
        Initialize the POM Simulator GUI application.

        This constructor sets up the main window of the application, including:
        - Window title and dimensions
        - Main layout with file system navigation panel
        - Tab structure for different simulation phases
        - Initial state for simulation threads

        The GUI is organized with a file browser on the left side and a tabbed
        interface on the right side containing all simulation functionality.

        Parameters:
            None

        Returns:
            None
        """
        super().__init__()
        self.setWindowTitle("POMSimulator")
        self.setGeometry(300, 200, 2048, 1152)

        # Set application logo
        logo_path = os.path.join(script_dir, "docs", ".img", "pomsimulator_logo.png")
        if os.path.exists(logo_path):
            self.setWindowIcon(QIcon(logo_path))
        self.pom_type = "IPA"
        # Initialize theme
        self.dark_theme = False
        self.font_size = 10.0

        # Initialize shared data dictionary for synchronization between tabs
        self.shared_data = {}

        # Initialize UI elements dictionary for tracking elements that share data
        self.ui_elements = {}

        # Create central widget with horizontal layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        self.create_mode_toolbar()
        self.create_menu()

        # Add file system dock on left
        self.create_file_dock()

        # Create main tab widget on right
        self.tabs = QTabWidget()
        main_layout.addWidget(self.tabs, 3)  # 3/4 of width

        # Create tabs

        self.refresh_tab_indicators()
        self.create_presimulation_tab()
        self.create_simulation_tab()
        self.create_scaling_tab()
        self.create_speciation_tab()
        self.create_clustering_tab()
        self.create_plotting_tab()
        # self.create_visualization_tab()

        self.tab_status = {
            'presimulation': {
                'generate_mol': 'idle',
                'compute_iso': 'idle'
            },
            'simulation': {
                'compute_lgkf': 'idle',
            },
            'scaling': 'idle',
            'speciation': {
                'spec_diag': 'idle',
                'phase_diag': 'idle'
            },
            'clustering': {
                'clustering': 'idle',
                'selection': 'idle',
                'filtering': 'idle'
            },
            'plotting': {
                'plot_spec': 'idle',
                'plot_phase': 'idle'
            }
        }

        self.refresh_tab_indicators()
        self.add_shared_console(main_layout)

        # Apply initial theme and font settings
        self.apply_theme()
        self.setup_shortcuts()

    def get_or_create_widget(self, attr_name, widget_class, *args, **kwargs):
        """
        Get an existing widget or create a new one if it doesn't exist.

        This method implements a widget registry system to prevent duplicate
        widget creation across different tab methods. It ensures that each
        widget attribute is created only once and reused across tabs.

        Args:
            attr_name (str): The attribute name for the widget (e.g., 'POM_system')
            widget_class: The PyQt5 widget class to instantiate (e.g., QLineEdit)
            *args: Positional arguments to pass to the widget constructor
            **kwargs: Keyword arguments to pass to the widget constructor

        Returns:
            The widget instance (either existing or newly created)

        Example:
            # Instead of: self.POM_system = QLineEdit()
            # Use: self.POM_system = self.get_or_create_widget('POM_system', QLineEdit)
        """
        if hasattr(self, attr_name):
            # Widget already exists, return the existing instance
            return getattr(self, attr_name)
        else:
            # Create new widget and store it
            widget = widget_class(*args, **kwargs)
            setattr(self, attr_name, widget)
            # Also track in ui_elements for potential future use
            self.ui_elements[attr_name] = widget
            return widget

    def register_shared_widget(self, attr_name, widget):
        """
        Register an existing widget in the registry system.

        This method allows manual registration of widgets that were created
        outside of get_or_create_widget() but should be tracked to prevent
        duplicate creation.

        Args:
            attr_name (str): The attribute name for the widget
            widget: The widget instance to register

        Returns:
            The registered widget instance
        """
        if not hasattr(self, attr_name):
            setattr(self, attr_name, widget)
            self.ui_elements[attr_name] = widget
        return getattr(self, attr_name)

    """Gui functions"""

    def create_menu(self):
        """
        Create the application menu bar with various options.

        This method sets up the menu bar for the application with File, Tools, and Help menus.
        It adds actions for opening files, accessing tools like the molecule visualizer,
        and displaying help information.

        Returns:
            None: This method modifies the GUI by adding a menu bar but doesn't return a value.
        """
        menubar = self.menuBar()

        # File menu
        file_menu = menubar.addMenu('File')

        # Add Open Configuration action
        open_action = QAction('Import Configuration File', self)
        open_action.setShortcut('Ctrl+O')
        open_action.setStatusTip("Open a configuration file")
        open_action.triggered.connect(self.load_config_file)
        file_menu.addAction(open_action)

        # Add Save Configuration action
        save_config_action = QAction("Export Configuration File...", self)
        save_config_action.setShortcut("Ctrl+S")
        save_config_action.setStatusTip("Save current configuration to a file")
        save_config_action.triggered.connect(self.save_config_file)
        file_menu.addAction(save_config_action)

        # Add separator and Exit action
        file_menu.addSeparator()
        exit_action = QAction("Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.setStatusTip("Exit the application")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        # Tools menu
        tools_menu = menubar.addMenu('Tools')

        # Add settings to Tools menu
        settings_action = QAction('Preferences', self)
        settings_action.setShortcut('Ctrl+P')
        settings_action.setToolTip("Open application preferences")
        settings_action.triggered.connect(self.open_settings)
        tools_menu.addAction(settings_action)

        # Add theme toggle to Tools menu
        theme_action = QAction('Toggle Theme', self)
        theme_action.setShortcut('Ctrl+T')
        theme_action.setToolTip("Switch between light and dark themes")
        theme_action.triggered.connect(self.toggle_theme)
        tools_menu.addAction(theme_action)

        mode_action = QAction('Toggle Mode', self)
        mode_action.setShortcut('Ctrl+M')
        mode_action.setToolTip("Switch between simulation and visualization modes")
        mode_action.triggered.connect(self.toggle_mode_with_shortcut)
        tools_menu.addAction(mode_action)

        refresh_gui = QAction('Refresh GUI', self)
        refresh_gui.setShortcut('Ctrl+F5')
        refresh_gui.setToolTip("Refresh the GUI to default parameters")
        refresh_gui.triggered.connect(self.refresh_interface)
        tools_menu.addAction(refresh_gui)

        # Help menu
        help_menu = menubar.addMenu('Help')

        docs_action = QAction('Documentation', self)
        docs_action.setShortcut('F2')
        docs_action.triggered.connect(self.open_documentation)
        help_menu.addAction(docs_action)

        about_action = QAction('About', self)
        about_action.setShortcut('Ctrl+A')
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)

        keybinding_action = QAction('Key Bindings', self)
        keybinding_action.setShortcut('Ctrl+F1')
        keybinding_action.triggered.connect(self.show_keybindings)
        help_menu.addAction(keybinding_action)

    def apply_theme(self):
        """
        Apply the current theme (light or dark) to the application.

        This method applies either the light or dark theme based on the current
        value of self.dark_theme. It ensures consistent styling between themes,
        only changing colors while maintaining the same visual structure.

        The light theme uses the default Qt styling with some minor customizations,
        while the dark theme uses a custom dark color palette that mirrors the
        structure of the light theme.

        Returns:
            None: This method modifies the application's appearance but doesn't return a value.
        """
        app = QApplication.instance()

        # Common stylesheet for both themes - tab sizing and font adjustments
        common_stylesheet = """
            QTabWidget::pane {
                border: 1px solid #c0c0c0;
                padding: 5px;
            }

            QTabBar::tab {
                font-size: 12pt;
                font-weight: bold;
                min-width: 200px;
                min-height: 30px;
                padding: 8px 12px;
                margin-right: 2px;
            }

            /* Style for nested tabs (smaller than parent tabs) */
            QTabWidget QTabWidget QTabBar::tab {
                font-size: 10pt;
                min-width: 250px;
                padding: 6px 10px;
            }

            /* Style for deeply nested tabs (even smaller) */
            QTabWidget QTabWidget QTabWidget QTabBar::tab {
                font-size: 8pt;
                min-width: 150px;
                padding: 4px 8px;
            }
        """

        if self.dark_theme:
            # Dark theme
            app.setStyle('Fusion')  # Use Fusion style for consistent appearance

            # Create a dark palette
            dark_palette = QPalette()

            # Set color groups
            dark_palette.setColor(QPalette.Window, QColor(53, 53, 53))
            dark_palette.setColor(QPalette.WindowText, QColor(255, 255, 255))
            dark_palette.setColor(QPalette.Base, QColor(35, 35, 35))
            dark_palette.setColor(QPalette.AlternateBase, QColor(45, 45, 45))
            dark_palette.setColor(QPalette.ToolTipBase, QColor(255, 255, 255))
            dark_palette.setColor(QPalette.ToolTipText, QColor(255, 255, 255))
            dark_palette.setColor(QPalette.Text, QColor(255, 255, 255))
            dark_palette.setColor(QPalette.Button, QColor(53, 53, 53))
            dark_palette.setColor(QPalette.ButtonText, QColor(255, 255, 255))
            dark_palette.setColor(QPalette.BrightText, QColor(255, 0, 0))
            dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
            dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
            dark_palette.setColor(QPalette.HighlightedText, QColor(255, 255, 255))

            # Set disabled colors
            dark_palette.setColor(QPalette.Disabled, QPalette.WindowText, QColor(127, 127, 127))
            dark_palette.setColor(QPalette.Disabled, QPalette.Text, QColor(127, 127, 127))
            dark_palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(127, 127, 127))
            dark_palette.setColor(QPalette.Disabled, QPalette.Highlight, QColor(80, 80, 80))
            dark_palette.setColor(QPalette.Disabled, QPalette.HighlightedText, QColor(127, 127, 127))

            # Apply the palette
            app.setPalette(dark_palette)

            # Dark theme specific tab styling
            dark_tab_style = """
                QTabBar::tab {
                    background-color: #3a3a3a;
                    color: #ffffff;
                    border: 1px solid #555555;
                    border-bottom: none;
                    border-top-left-radius: 4px;
                    border-top-right-radius: 4px;
                }

                QTabBar::tab:selected {
                    background-color: #505050;
                    border-bottom: none;
                }

                QTabBar::tab:hover:!selected {
                    background-color: #454545;
                }
            """

            # Apply combined stylesheet
            app.setStyleSheet(common_stylesheet + dark_tab_style)

        else:
            # Light theme - reset to default
            app.setStyle('Fusion')  # Use Fusion style for consistent appearance
            app.setPalette(app.style().standardPalette())  # Reset to default palette

            # Light theme specific tab styling
            light_tab_style = """
                QTabBar::tab {
                    background-color: #f0f0f0;
                    color: #000000;
                    border: 1px solid #c0c0c0;
                    border-bottom: none;
                    border-top-left-radius: 4px;
                    border-top-right-radius: 4px;
                }

                QTabBar::tab:selected {
                    background-color: #ffffff;
                    border-bottom: none;
                }

                QTabBar::tab:hover:!selected {
                    background-color: #e0e0e0;
                }
            """

            # Apply combined stylesheet
            app.setStyleSheet(common_stylesheet + light_tab_style)

        self.update_console_theme()  # Update the console theme to match the new theme

    def toggle_theme(self):
        """
        Toggle between light and dark themes for the application.

        This method switches the application's stylesheet between light and dark themes.
        It updates the self.dark_theme flag to track the current theme state and applies
        the appropriate stylesheet to the entire application.

        Returns:
            None: This method modifies the application's appearance but doesn't return a value.
        """
        # Toggle the theme state
        self.dark_theme = not self.dark_theme

        # Apply the theme without updating the theme selector
        self._updating_theme = True  # Set a flag to prevent recursion
        self.apply_theme()
        self.update_console_theme()
        self.update_constants_display()
        self._updating_theme = False

        # Update theme selection in settings dialog if it's open
        if hasattr(self, 'theme_selector') and self.theme_selector is not None:
            # Temporarily disconnect the signal to prevent recursion
            try:
                self.theme_selector.blockSignals(True)
                self.theme_selector.setCurrentIndex(1 if self.dark_theme else 0)
            finally:
                self.theme_selector.blockSignals(False)

    def show_about(self):
        """
        Show the About dialog with information about POMSimulator.
        """
        about_text = """
        <h2>POMSimulator 2.0</h2>
        <p>A comprehensive simulation tool for aqueous polyoxometalate speciation.</p>
        <p>Authors: ENRIC PETRUS, MIREIA SEGADO-CENTELLAS, CARLES BO.</p>
        <p>Developed by the ENRIC PETRUS, JORDI BUILS, DIEGO GARAY-RUIZ.</p>
        <p>Version: 2.0</p>
        <p>POMSimulator is distributed under the GNU AFFERO GENERAL PUBLIC LICENSE</p>
        """
        msg_box = QMessageBox(self)
        msg_box.setWindowTitle("About POMSimulator")
        msg_box.setText(about_text)
        msg_box.setIcon(QMessageBox.Information)
        msg_box.setStandardButtons(QMessageBox.Ok)
        msg_box.resize(600, 300)
        msg_box.exec_()

    def show_keybindings(self):
        """
        Show a dialog with all available keyboard shortcuts and their descriptions.
        """
        try:
            # Determine colors based on current theme
            if self.dark_theme:
                # Dark theme colors
                header_bg = "#404040"
                subheader_bg = "#505050"
                row_bg = "#353535"
                alt_row_bg = "#454545"
                text_color = "#ffffff"
                border_color = "#666666"
            else:
                # Light theme colors
                header_bg = "#d0d0d0"
                subheader_bg = "#f0f0f0"
                row_bg = "#ffffff"
                alt_row_bg = "#f9f9f9"
                text_color = "#000000"
                border_color = "#c0c0c0"

            # Create the keybindings information organized by menu
            keybindings_text = f"""
            <h2 style="color: {text_color};">Keyboard Shortcuts</h2>
            <table border="1" cellpadding="8" cellspacing="0" style="border-collapse: collapse; width: 100%; border-color: {border_color};">

            <!-- File Menu Section -->
            <tr style="background-color: {header_bg};">
                <th colspan="3" style="text-align: left; font-weight: bold; font-size: 12pt; color: {text_color};">File Menu</th>
            </tr>
            <tr style="background-color: {subheader_bg};">
                <th style="text-align: left; font-weight: bold; color: {text_color};">Shortcut</th>
                <th style="text-align: left; font-weight: bold; color: {text_color};">Action</th>
                <th style="text-align: left; font-weight: bold; color: {text_color};">Description</th>
            </tr>
            <tr style="background-color: {row_bg};">
                <td style="color: {text_color};"><b>Ctrl+O</b></td>
                <td style="color: {text_color};">Import Configuration</td>
                <td style="color: {text_color};">Load parameters from configuration file</td>
            </tr>
            <tr style="background-color: {alt_row_bg};">
                <td style="color: {text_color};"><b>Ctrl+S</b></td>
                <td style="color: {text_color};">Export Configuration</td>
                <td style="color: {text_color};">Save current parameters to configuration file</td>
            </tr>
            <tr style="background-color: {row_bg};">
                <td style="color: {text_color};"><b>Ctrl+Q</b></td>
                <td style="color: {text_color};">Exit Application</td>
                <td style="color: {text_color};">Close POMSimulator</td>
            </tr>

            <!-- Tools Menu Section -->
            <tr style="background-color: {header_bg};">
                <th colspan="3" style="text-align: left; font-weight: bold; font-size: 12pt; color: {text_color};">Tools Menu</th>
            </tr>
            <tr style="background-color: {subheader_bg};">
                <th style="text-align: left; font-weight: bold; color: {text_color};">Shortcut</th>
                <th style="text-align: left; font-weight: bold; color: {text_color};">Action</th>
                <th style="text-align: left; font-weight: bold; color: {text_color};">Description</th>
            </tr>
            <tr style="background-color: {row_bg};">
                <td style="color: {text_color};"><b>Ctrl+P</b></td>
                <td style="color: {text_color};">Preferences</td>
                <td style="color: {text_color};">Open application preferences dialog</td>
            </tr>
            <tr style="background-color: {alt_row_bg};">
                <td style="color: {text_color};"><b>Ctrl+T</b></td>
                <td style="color: {text_color};">Toggle Theme</td>
                <td style="color: {text_color};">Switch between light and dark themes</td>
            </tr>
            <tr style="background-color: {row_bg};">
                <td style="color: {text_color};"><b>Ctrl++</b></td>
                <td style="color: {text_color};">Increase Text Size</td>
                <td style="color: {text_color};">Make text larger throughout the application</td>
            </tr>
            <tr style="background-color: {alt_row_bg};">
                <td style="color: {text_color};"><b>Ctrl+-</b></td>
                <td style="color: {text_color};">Decrease Text Size</td>
                <td style="color: {text_color};">Make text smaller throughout the application</td>
            </tr>

            <!-- Help Menu Section -->
            <tr style="background-color: {header_bg};">
                <th colspan="3" style="text-align: left; font-weight: bold; font-size: 12pt; color: {text_color};">Help Menu</th>
            </tr>
            <tr style="background-color: {subheader_bg};">
                <th style="text-align: left; font-weight: bold; color: {text_color};">Shortcut</th>
                <th style="text-align: left; font-weight: bold; color: {text_color};">Action</th>
                <th style="text-align: left; font-weight: bold; color: {text_color};">Description</th>
            </tr>
            <tr style="background-color: {row_bg};">
                <td style="color: {text_color};"><b>F2</b></td>
                <td style="color: {text_color};">Documentation</td>
                <td style="color: {text_color};">Open POMSimulator documentation in browser</td>
            </tr>
            <tr style="background-color: {alt_row_bg};">
                <td style="color: {text_color};"><b>Ctrl+A</b></td>
                <td style="color: {text_color};">About</td>
                <td style="color: {text_color};">Show information about POMSimulator</td>
            </tr>
            <tr style="background-color: {row_bg};">
                <td style="color: {text_color};"><b>Ctrl+F1</b></td>
                <td style="color: {text_color};">Key Bindings</td>
                <td style="color: {text_color};">Show this keyboard shortcuts dialog</td>
            </tr>

            <!-- Interface Navigation Section -->
            <tr style="background-color: {header_bg};">
                <th colspan="3" style="text-align: left; font-weight: bold; font-size: 12pt; color: {text_color};">Interface Navigation</th>
            </tr>
            <tr style="background-color: {subheader_bg};">
                <th style="text-align: left; font-weight: bold; color: {text_color};">Shortcut</th>
                <th style="text-align: left; font-weight: bold; color: {text_color};">Action</th>
                <th style="text-align: left; font-weight: bold; color: {text_color};">Description</th>
            </tr>
            <tr style="background-color: {row_bg};">
                <td style="color: {text_color};"><b>Ctrl+M</b></td>
                <td style="color: {text_color};">Toggle Mode</td>
                <td style="color: {text_color};">Switch between IPA and HPA modes</td>
            </tr>
            <tr style="background-color: {alt_row_bg};">
                <td style="color: {text_color};"><b>Ctrl+F5</b></td>
                <td style="color: {text_color};">Refresh Interface</td>
                <td style="color: {text_color};">Clear all parameters and refresh the interface</td>
            </tr>

            </table>
            <br>
            <p style="color: {text_color};"><i>Note: Some shortcuts may only be available in specific contexts or modes.</i></p>
                    """

            # Create and show the message box
            msg_box = QMessageBox(self)
            msg_box.setWindowTitle("Keyboard Shortcuts - POMSimulator")
            msg_box.setTextFormat(Qt.RichText)
            msg_box.setText(keybindings_text)
            msg_box.setIcon(QMessageBox.Information)
            msg_box.setStandardButtons(QMessageBox.Ok)

            # Make the dialog larger to accommodate the table
            msg_box.resize(700, 700)

            # Show the dialog
            msg_box.exec_()

        except Exception as e:
            print(f"Error showing keybindings dialog: {e}")
            # Fallback to simple message if rich text fails
            QMessageBox.information(
                self,
                "Keyboard Shortcuts",
                "FILE MENU:\n"
                "Ctrl+O: Import configuration\n"
                "Ctrl+S: Export configuration\n"
                "Ctrl+Q: Exit application\n\n"
                "TOOLS MENU:\n"
                "Ctrl+P: Preferences\n"
                "Ctrl+T: Toggle theme\n"
                "Ctrl++: Increase text size\n"
                "Ctrl+-: Decrease text size\n\n"
                "HELP MENU:\n"
                "F2: Documentation\n"
                "Ctrl+A: About\n"
                "Ctrl+F1: Key bindings\n\n"
                "INTERFACE:\n"
                "Ctrl+M: Toggle IPA/HPA modes\n"
                "Ctrl+F5: Refresh interface"
            )

    def open_documentation(self):
        """
        Open the POMSimulator documentation in the default web browser.

        This method uses QDesktopServices to open the documentation URL in the user's
        default web browser. If the operation fails, it displays an error message.

        Returns:
            None: This method opens a web browser but doesn't return a value.
        """
        from PyQt5.QtCore import QUrl

        # URL to documentation
        doc_url = "https://pomsimulator.readthedocs.io/en/latest/"

        try:
            QDesktopServices.openUrl(QUrl(doc_url))
        except Exception as e:
            QMessageBox.warning(self, "Error Opening Documentation",
                                f"Could not open documentation: {str(e)}")

    def open_settings(self):
        """
        Open a simplified settings dialog with appearance options.

        This dialog allows users to:
        - Select between light and dark themes
        - Adjust text size with buttons (with keyboard shortcuts)
        """
        # Create dialog
        dialog = QDialog(self)
        dialog.setWindowTitle("Appearance Settings")
        dialog.setMinimumWidth(400)

        # Create layout
        layout = QVBoxLayout(dialog)

        # Theme selection group
        theme_group = QGroupBox("Theme")
        theme_layout = QHBoxLayout()

        # Create theme radio buttons
        light_theme = QRadioButton("Light")
        dark_theme = QRadioButton("Dark")

        # Set the current theme
        if self.dark_theme:
            dark_theme.setChecked(True)
        else:
            light_theme.setChecked(True)

        theme_layout.addWidget(light_theme)
        theme_layout.addWidget(dark_theme)
        theme_group.setLayout(theme_layout)
        layout.addWidget(theme_group)

        # Text size group
        text_group = QGroupBox("Text Size")
        text_layout = QHBoxLayout()

        # Current size label
        current_size_label = QLabel("Current Size:")
        self.text_size_value = QLabel(f"{self.font_size:.1f}pt")

        # Create buttons for text size adjustment
        decrease_btn = QPushButton("Smaller (Ctrl+-)")
        decrease_btn.setShortcut("Ctrl+-")
        decrease_btn.clicked.connect(lambda: self.adjust_text_size(-0.5, self.text_size_value))

        increase_btn = QPushButton("Larger (Ctrl++)")
        increase_btn.setShortcut("Ctrl++")
        increase_btn.clicked.connect(lambda: self.adjust_text_size(0.5, self.text_size_value))

        text_layout.addWidget(current_size_label)
        text_layout.addWidget(self.text_size_value)
        text_layout.addWidget(decrease_btn)
        text_layout.addWidget(increase_btn)
        text_group.setLayout(text_layout)
        layout.addWidget(text_group)

        # Add buttons at the bottom
        button_layout = QHBoxLayout()
        apply_btn = QPushButton("Apply")
        apply_btn.clicked.connect(lambda: self.apply_settings(
            dark_theme.isChecked(),
            dialog
        ))
        button_layout.addWidget(apply_btn)

        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(dialog.reject)
        button_layout.addWidget(cancel_btn)

        layout.addLayout(button_layout)

        # Show dialog
        dialog.exec_()

    def adjust_text_size(self, delta, label=None):
        """
        Adjust the application's text size.

        Parameters:
            delta (float): Amount to change the font size by
            label (QLabel, optional): Label to update with new size
        """
        # Initialize font_size attribute if it doesn't exist
        if not hasattr(self, 'font_size'):
            self.font_size = 10.0

        # Adjust the font size
        new_size = self.font_size + delta

        # Limit the size range
        if 6.0 <= new_size <= 16.0:
            self.font_size = new_size

            # Update the label if provided
            if label:
                label.setText(f"{self.font_size:.1f}pt")

            # Apply the new font size to the application
            self.apply_font_size()

    def apply_font_size(self):
        """
        Apply the current font size to the application.
        """
        app = QApplication.instance()

        # Get the current stylesheet
        current_style = app.styleSheet()

        # Remove any existing font-size declarations
        style_lines = current_style.split('\n')
        filtered_lines = []
        for line in style_lines:
            if 'font-size:' not in line:
                filtered_lines.append(line)

        # Add the new font size
        font_style = f"""
        * {{
            font-size: {self.font_size:.1f}pt;
        }}
        """

        # Combine and apply
        new_style = '\n'.join(filtered_lines) + font_style
        app.setStyleSheet(new_style)

    def apply_settings(self, dark_mode, dialog):
        """
        Apply the selected settings and close the dialog.

        Parameters:
            dark_mode (bool): Whether dark mode is selected
            dialog (QDialog): The settings dialog to close
        """
        # Apply theme if changed
        if self.dark_theme != dark_mode:
            self.dark_theme = dark_mode
            self.apply_theme()

        # Apply font size
        self.apply_font_size()

        # Close dialog
        dialog.accept()

    def setup_shortcuts(self):
        """
        Set up keyboard shortcuts for the application.
        """
        # Text size shortcuts
        increase_shortcut = QShortcut(QKeySequence("Ctrl++"), self)
        increase_shortcut.activated.connect(lambda: self.adjust_text_size(0.5))

        decrease_shortcut = QShortcut(QKeySequence("Ctrl+-"), self)
        decrease_shortcut.activated.connect(lambda: self.adjust_text_size(-0.5))

    def create_file_dock(self):
        """
        Create a file navigation dock widget for the main window.

        This method sets up a dockable file browser panel on the left side of the
        application window. It creates a file system model rooted at the pomsimulator
        directory and displays it in a tree view, allowing users to navigate the
        project's file structure.

        The dock widget can be moved between the left and right sides of the main window.
        The tree view shows only file names for cleaner navigation. Double-clicking on files
        will open them with the appropriate application based on their file type:
        - Text files will be opened in a text viewer within the application
        - Image files will be displayed in an image viewer
        - Other file types will be opened with the system's default application

        Parameters:
            None

        Returns:
            None
        """
        dock = QDockWidget("File System", self)
        dock.setAllowedAreas(Qt.LeftDockWidgetArea | Qt.RightDockWidgetArea)

        # Determine the correct root path for both development and compiled environments
        def get_project_root():
            """Get the project root directory, handling both development and PyInstaller environments."""
            if getattr(sys, 'frozen', False):
                # Running as compiled executable
                # Get the directory where the executable is located
                executable_dir = os.path.dirname(sys.executable)

                # Look for pomsimulator directory structure
                # Check if we're in the pomsimulator directory already
                if os.path.basename(executable_dir) == 'pomsimulator':
                    return executable_dir

                # Check if pomsimulator is a subdirectory
                pomsim_path = os.path.join(executable_dir, 'pomsimulator')
                if os.path.exists(pomsim_path):
                    return pomsim_path

                # If not found, use the executable directory as fallback
                return executable_dir
            else:
                # Running in development mode
                return script_dir

        root_path = get_project_root()
        parent_path = root_path
        if "pomsimulator" not in root_path:
            parent_path = os.path.dirname(root_path)

        # Create file system model
        self.model = QFileSystemModel()
        self.model.setRootPath(os.path.join(parent_path))

        # Create tree view
        self.tree = QTreeView()
        self.tree.setModel(self.model)
        self.tree.setRootIndex(self.model.index(os.path.join(parent_path)))

        # Hide all columns except the first one (name)
        for i in range(1, self.model.columnCount()):
            self.tree.hideColumn(i)

        # Set appropriate width for the name column
        self.tree.setColumnWidth(0, 450)

        # Connect double-click event to file opening function
        self.tree.doubleClicked.connect(self.open_file_from_browser)

        dock.setWidget(self.tree)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock)

    def open_file_from_browser(self, index):
        """
        Open a file when it's double-clicked in the file browser.

        This method is triggered when a user double-clicks on an item in the file browser.
        If the clicked item is a file (not a directory), it determines the file type and
        opens it appropriately:
        - Text files (.txt, .py, .csv, etc.) are opened in a text viewer dialog
        - Image files (.png, .jpg, etc.) are opened in an image viewer dialog
        - Molecule files (.mol) are opened in the molecule visualizer
        - Other file types are opened using the system's default application

        Parameters:
            index (QModelIndex): The index of the item that was double-clicked in the
                            file browser tree view.

        Returns:
            None: This method opens a file or displays a dialog but doesn't return a value.
        """
        # Get the file path from the model index
        file_path = self.model.filePath(index)

        # Check if it's a file (not a directory)
        if os.path.isfile(file_path):
            file_ext = os.path.splitext(file_path)[1].lower()

            # Handle text files
            if file_ext in ['.txt', '.py', '.csv', '.json', '.xml', '.md', '.log', '.ini', '.pomsim', '', '.out']:
                self.open_text_file(file_path)

            # Handle image files
            elif file_ext in ['.png', '.jpg', '.jpeg', '.gif', '.bmp', '.tiff', '.svg']:
                self.open_image_file(file_path)

            # Handle molecule files
            elif file_ext == '.mol':
                self.open_mol_file(file_path)

            elif file_ext == '.ipynb':
                QMessageBox.warning(self, "Error Opening File", "No viewer available for this file type.")

            # Handle other file types with system default application
            else:
                try:
                    # Use the appropriate method based on the platform
                    if sys.platform.startswith('darwin'):  # macOS
                        subprocess.call(('open', file_path))
                    elif os.name == 'nt':  # Windows
                        os.startfile(file_path)
                    elif os.name == 'posix':  # Linux
                        subprocess.call(('xdg-open', file_path))
                except Exception as e:
                    QMessageBox.warning(self, "Error Opening File",
                                        f"Could not open file: {str(e)}")

    def open_text_file(self, file_path):
        """
        Open a text file in a viewer dialog.

        This method creates a dialog window to display the contents of a text file.
        The dialog includes the file name in its title and displays the text content
        in a read-only text edit widget with a monospaced font for better readability.

        Parameters:
            file_path (str): The full path to the text file to be opened.

        Returns:
            None: This method displays a dialog but doesn't return a value.
        """
        try:
            with open(file_path, 'r', encoding='utf-8') as file:
                content = file.read()

            # Create dialog
            dialog = QDialog(self)
            dialog.setWindowTitle(f"Text Viewer - {os.path.basename(file_path)}")
            dialog.resize(800, 800)

            # Create layout
            layout = QVBoxLayout(dialog)

            # Create text edit
            text_edit = QTextEdit()
            text_edit.setReadOnly(False)
            if file_path.endswith('.out') or file_path.endswith('.log'):
                text_edit.setReadOnly(True)
            text_edit.setPlainText(content)

            # Set monospace font
            monospace_font = QFont("Monospace")
            monospace_font.setStyleHint(QFont.Monospace)
            monospace_font.setPointSize(10)
            text_edit.setFont(monospace_font)

            layout.addWidget(text_edit)

            # Add button layout
            button_layout = QHBoxLayout()

            # Add save button
            save_btn = QPushButton("Save")
            save_btn.clicked.connect(lambda: self.save_text_file(file_path, text_edit.toPlainText()))
            button_layout.addWidget(save_btn)

            # Add close button
            close_btn = QPushButton("Close")
            close_btn.clicked.connect(dialog.close)
            button_layout.addWidget(close_btn)

            layout.addLayout(button_layout)

            # Add keyboard shortcut for saving
            save_shortcut = QShortcut(QKeySequence("Ctrl+S"), dialog)
            save_shortcut.activated.connect(lambda: self.save_text_file(file_path, text_edit.toPlainText()))

            dialog.setLayout(layout)
            dialog.exec_()

        except Exception as e:
            QMessageBox.warning(self, "Error Opening Text File",
                                f"Could not open text file: {str(e)}")

    def save_text_file(self, file_path, content):
        """
        Save the modified text content back to the file.

        Parameters:
            file_path (str): The full path to the text file to be saved.
            content (str): The text content to save.

        Returns:
            None
        """
        try:
            with open(file_path, 'w', encoding='utf-8') as file:
                file.write(content)
            QMessageBox.information(self, "File Saved", f"File saved successfully: {os.path.basename(file_path)}")
        except Exception as e:
            QMessageBox.warning(self, "Error Saving File", f"Could not save file: {str(e)}")

    def open_image_file(self, file_path):
        """
        Open an image file in a viewer dialog.

        This method creates a dialog window to display an image file. The dialog
        includes the file name in its title and shows the image in a label widget
        that's sized appropriately for viewing.

        Parameters:
            file_path (str): The full path to the image file to be opened.

        Returns:
            None: This method displays a dialog but doesn't return a value.
        """
        try:
            # Create dialog
            dialog = QDialog(self)
            dialog.setWindowTitle(f"Image Viewer - {os.path.basename(file_path)}")
            dialog.resize(800, 800)

            # Create layout
            layout = QVBoxLayout(dialog)

            # Create scroll area for the image
            scroll_area = QScrollArea()
            scroll_area.setWidgetResizable(True)

            # Create image label
            image_label = QLabel()
            pixmap = QPixmap(file_path)
            image_label.setPixmap(pixmap)
            image_label.setAlignment(Qt.AlignCenter)

            # Add image to scroll area
            scroll_area.setWidget(image_label)
            layout.addWidget(scroll_area)

            # Add close button
            close_btn = QPushButton("Close")
            close_btn.clicked.connect(dialog.close)
            layout.addWidget(close_btn)

            dialog.setLayout(layout)
            dialog.exec_()

        except Exception as e:
            QMessageBox.warning(self, "Error Opening Image File",
                                f"Could not open image file: {str(e)}")

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
            from utilities.ase_visualizer import MoleculeVisualizer

            # Create the visualizer
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

    def save_config_file(self, file_path=None):
        """
        Save the current GUI parameters to a configuration file..

        This method creates a dialog window to select a file to save the configuration
        settings. The dialog includes the file name in its title and allows users to
        choose a location and name for the configuration file.This method collects all
        relevant parameters from the GUI and writes them to a configuration file in the
        *.pomsim format. It handles both IPA (Isopolyanions) and HPA (Heteropolyanions)
        configurations differently based on the POM system name.

        Parameters:
            file_path (str, optional): The full path to the configuration file to be saved.
                                    If not provided, a default file name will be used.

        Returns:
            bool: True if the configuration was saved successfully, False otherwise.
        """
        try:
            # If no file path is provided, open a file dialog
            if not file_path:
                file_path, _ = QFileDialog.getSaveFileName(self, "Save Configuration File", os.getcwd(),
                                                           "Configuration Files (*.pomsim);;All Files (*)",
                                                           options=QFileDialog.DontUseNativeDialog)
            if not file_path:
                return False

            # Add .pomsim extension if not already present
            if not file_path.lower().endswith('.pomsim'):
                file_path += '.pomsim'
        except Exception as e:
            self.console.append(f"Error: Failed to save configuration file: {str(e)}")
            QMessageBox.warning(self, "Error Opening File", f"Failed to save configuration file: {str(e)}")
            return False

        try:
            config_dict = self.collect_config()

            config = ConfigParser()
            config.optionxform = str
            for section, values in config_dict.items():
                config.add_section(section)
                for key, value in values.items():
                    config.set(section, key, str(value))

            with open(file_path, 'w') as configfile:
                config.write(configfile)

        except Exception as e:
            self.update_console("Error: Failed to save configuration file.", 'presimulation')
            QMessageBox.warning(self, "Error Saving Configuration",
                                f"Failed to save configuration file: {str(e)}")
            return False

    def load_config_file(self):
        """
        Load configuration from a file and apply settings to the GUI.

        Returns:
            bool: True if configuration was loaded successfully, False otherwise.
        """
        try:
            config_file_path, _ = QFileDialog.getOpenFileName(
                self,
                "Open Configuration File",
                os.path.join(script_dir, "inputs"),
                "Configuration Files (*.pomsim);;All Files (*)",
                options=QFileDialog.DontUseNativeDialog
            )

            if not config_file_path:
                return False

            # Define the set_widget_value method to handle different widget types
            import types

            def set_widget_value(self, widget, value):
                """
                Set the value of a widget based on its type.

                Args:
                    widget: The widget to set the value for
                    value: The value to set
                """
                try:
                    if isinstance(widget, QLineEdit):
                        widget.setText(str(value))
                    elif isinstance(widget, QSpinBox):
                        try:
                            # Ensure we're passing an integer to QSpinBox
                            int_value = int(value)
                            widget.setValue(int_value)
                        except (ValueError, TypeError) as e:
                            self.console.append(f"Warning: Could not convert '{value}' to an integer: {str(e)}")
                    elif isinstance(widget, CustomDoubleSpinBox):
                        try:
                            # Convert to float for QDoubleSpinBox
                            float_value = float(value)
                            widget.setValue(float_value)
                        except (ValueError, TypeError) as e:
                            self.console.append(f"Warning: Could not convert '{value}' to a float: {str(e)}")
                    elif isinstance(widget, QComboBox):
                        index = widget.findText(str(value))
                        if index >= 0:
                            widget.setCurrentIndex(index)
                        else:
                            # Try case-insensitive match
                            for i in range(widget.count()):
                                if widget.itemText(i).lower() == str(value).lower():
                                    widget.setCurrentIndex(i)
                                    break
                    elif isinstance(widget, QCheckBox):
                        # Convert string to boolean
                        if isinstance(value, str):
                            widget.setChecked(value.lower() == "true")
                        else:
                            widget.setChecked(bool(value))
                except Exception as e:
                    self.update_console(f"Error updating {widget.__class__.__name__} for element {widget}: {str(e)}",
                                        'presimulation')

            # Attach the method to the class instance
            self.set_widget_value = types.MethodType(set_widget_value, self)

            config = ConfigParser()
            config.read(config_file_path)

            # Determine if the system is IPA or HPA based on the POM_system parameter
            # system_type = "ipa"  # Default to IPA
            if config.has_section('Preparation') and config.has_option('Preparation', 'POM_system'):
                pom_system = config.get('Preparation', 'POM_system')
                if "_" in pom_system:  # HPA systems typically have an underscore in their name
                    system_type = "hpa"
                else:  # IPA systems typically do not have an underscore in their name
                    system_type = "ipa"

            QApplication.processEvents()

            required_mode = "HPA" if system_type == "hpa" else "IPA"

            current_mode = getattr(self, 'global_mode', 'IPA')
            if hasattr(self, 'global_mode') and hasattr(self.global_mode, 'currentText'):
                current_mode = self.global_mode.currentText()

            # Switch to the appropriate tab based on system type
            if current_mode != required_mode:
                # Show mode change confirmation
                mode_name = "Heteropolyoxometalates" if required_mode == "HPA" else "Isopolyoxometalates"
                current_mode_name = "Heteropolyoxometalates" if current_mode == "HPA" else "Isopolyoxometalates"

                reply = QMessageBox.question(
                    self,
                    "Mode Change Required",
                    f"The configuration file requires {mode_name} mode, but you are currently in {current_mode_name} mode.\n\n"
                    f"Do you want to switch to {mode_name} mode and load the configuration?",
                    QMessageBox.Yes | QMessageBox.No,
                    QMessageBox.Yes
                )

                if reply == QMessageBox.Yes:
                    # Execute mode change using the new system
                    self.execute_mode_change(required_mode)
                    QApplication.processEvents()  # Process events to ensure UI updates
                else:
                    # User declined mode change
                    QMessageBox.information(
                        self,
                        "Configuration Load Cancelled",
                        "Configuration loading was cancelled due to mode mismatch."
                    )
                    return False

            # Define path fields that should be transformed from absolute to relative
            path_fields = {
                'Preparation': ['adf_inputs_dir', 'mol_folder', 'output_path'],
                'Speciation': ['npz_file', 'phase_dir', 'model_subset_file'],
                'Clustering': ['cluster_dir', 'npz_cluster_file', 'features_file'],
                'Visualization': ['labels_file']
            }
            # Map configuration sections to widget attributes
            # This dictionary maps config keys to widget attribute names
            widget_mapping = {
                'Preparation': {
                    'POM_system': 'POM_system',
                    'adf_inputs_dir': 'adf_inputs_dir',
                    'mol_folder': 'mol_folder',
                    'output_path': 'output_path'
                },
                'Isomorphism': {
                    'cores': 'iso_cores'
                },
                'Simulation': {
                    'cores': 'sim_cores',
                    'batch_size': 'sim_batch_size',
                    'sample_perc': 'sim_sample_perc',
                    'sample_type': 'sim_sample_type',
                    'use_isomorphism': 'use_isomorphism',
                    'energy_threshold': 'energy_threshold',
                    'proton_numb': 'proton_numb',
                    'reference_types': 'ref_types',
                    'ref_types': 'ref_types',
                    'I': 'i_s',
                    'C0': 'c0',
                    'temp': 'temp',
                    'min_pH': 'sim_min_ph',
                    'max_pH': 'sim_max_ph',
                    'step_pH': 'sim_step_ph',
                    'ref_compound': 'ref_compound',
                    'CM': 'sim_CM',
                    'CX': 'sim_CX'
                },
                'CRN': {
                    'Full_CRN': 'full_crn',
                    'Selected_model': 'selected_model',
                    'Plot_3D': 'plot_3d'
                },
                'Scaling': {
                    'scaling_mode': 'scaling_mode',
                    'experimental_set': 'exp_set'
                },
                'Speciation': {
                    'speciation_labels': 'spec_labels',
                    'min_pH': 'spec_min_ph',
                    'max_pH': 'spec_max_ph',
                    'step_pH': 'spec_step_ph',
                    'cores': 'spec_cores',
                    'batch_size': 'spec_batch_size',
                    'm_idx': 'm_idx',
                    'phase_dir': 'phase_dir',
                    'model_subset_file': 'model_subset_file',
                    'C': 'spec_C0',
                    'min_logC': 'phase_min_c',
                    'max_logC': 'phase_max_c',
                    'num_logC': 'phase_num_c',
                    'C_M': 'spec_CM',
                    'C_X': 'spec_CX',
                    'min_Ratio': 'min_metal_ratio',
                    'max_Ratio': 'max_metal_ratio',
                    'num_Ratio': 'num_metal_ratio'
                },
                'InternalConditions': {
                    'restrain_addition': 'rest_add',
                    'restrain_condensation': 'rest_cond',
                    'include_dimerization': 'include_dimerization',
                    'force_stoich': 'force_sto',
                    'adjust_protons_hydration': 'adj_prot'
                },
                'Clustering': {
                    'cluster_dir': 'cluster_dir',
                    'npz_cluster_file': 'npz_cluster_file',
                    'features_file': 'features_file',
                    'n_clusters': 'n_clusters',
                    'normalize_feats': 'normalize_feats',
                    'feats_list': 'feats_list',
                    'sel_groups': 'sel_groups'
                },
                'Visualization': {
                    'col_dict': 'clustering_col_dict.col_dict_combo',
                    'plot_list': 'plot_list',
                    'boxplot_list': 'boxplot_list'
                }
            }
            # First, handle the POM_system field specifically to ensure it's updated
            if config.has_section('Preparation') and config.has_option('Preparation', 'POM_system'):
                pom_system_value = config.get('Preparation', 'POM_system')
                # Try different attribute names that might hold the POM system field
                for attr_name in ['POM_system']:
                    if hasattr(self, attr_name):
                        widget = getattr(self, attr_name)
                        self.set_widget_value(widget, pom_system_value)
                        # self.update_console(f"Set {attr_name} to {pom_system_value}","Presimulation")
                        break
            # Load values into form fields
            for section in config.sections():
                if section in widget_mapping:
                    for option in config.options(section):
                        if option in widget_mapping[section]:
                            widget_name = widget_mapping[section][option]
                            if '.' in widget_name:
                                parts = widget_name.split('.')
                                if hasattr(self, parts[0]):
                                    parent_widget = getattr(self, parts[0])
                                    if hasattr(parent_widget, parts[1]):
                                        widget = getattr(parent_widget, parts[1])
                                        value = config.get(section, option)
                                        if section in path_fields and option in path_fields[section]:
                                            value = self.transform_absolute_to_relative_path(value)
                                        self.set_widget_value(widget, value)
                            else:
                                if hasattr(self, widget_name):
                                    widget = getattr(self, widget_name)
                                    value = config.get(section, option)

                                    # self.update_console(f"Set {widget_name} to {value}","Presimulation")
                                    if section in path_fields and option in path_fields[section]:
                                        value = self.transform_absolute_to_relative_path(value)
                                    self.set_widget_value(widget, value)

            # Special handling for reference compounds in HPA mode
            if system_type == "hpa" and config.has_section('Simulation') and config.has_option('Simulation',
                                                                                               'ref_compound'):
                ref_compounds = config.get('Simulation', 'ref_compound').split(',')
                if len(ref_compounds) >= 1 and hasattr(self, 'sim_M_ref_compound'):
                    self.set_widget_value(self.sim_M_ref_compound, ref_compounds[0])
                if len(ref_compounds) >= 2 and hasattr(self, 'sim_X_ref_compound'):
                    self.set_widget_value(self.sim_X_ref_compound, ref_compounds[1])
            if not config.has_option('Visualization', 'plot_list') and hasattr(self, 'plot_list'):
                self.plot_list.setText("all")

            QMessageBox.information(self, "Configuration Loaded", "Configuration loaded successfully.")
            # self.update_console(f"Configuration loaded from {config_file_path}","Presimulation")
            return True

        except Exception as e:
            import traceback
            error_msg = f"Error loading configuration: {str(e)}\n{traceback.format_exc()}"
            self.update_console(error_msg, "Presimulation")
            QMessageBox.warning(self, "Error Loading Configuration", f"Failed to load configuration: {str(e)}")
            return False

    def collect_config(self):
        """
        Collect all current GUI parameter values into a nested dictionary.

        This is the single source of truth for reading the GUI state. Both
        ``save_config_file()`` and ``get_gui_params()`` delegate to this method
        so that adding or renaming a parameter only requires a change here.

        Returns:
            dict: Nested dictionary with section names as keys and parameter
                  dicts as values, matching the structure of the .pomsim config file.
        """
        system_name = self.POM_system.text() if hasattr(self, 'POM_system') else ""
        is_ipa = "_" not in system_name

        def get_shared_line_text(attribute_name, shared_key):
            """Return synchronized QLineEdit text, preferring shared_data over stale tab attributes."""
            shared_value = self.shared_data.get(shared_key, "") if hasattr(self, 'shared_data') else ""
            if shared_value not in (None, ""):
                return str(shared_value)
            widget = getattr(self, attribute_name, None)
            return widget.text() if isinstance(widget, QLineEdit) else ""

        config_dict = {
            "Preparation": {
                'POM_system': get_shared_line_text('POM_system', 'POM_system'),
                'adf_inputs_dir': get_shared_line_text('adf_inputs_dir', 'adf_inputs_dir'),
                'mol_folder': get_shared_line_text('mol_folder', 'mol_folder'),
                'output_path': get_shared_line_text('output_path', 'output_path'),
            },
            "Isomorphism": {
                'cores': str(self.iso_cores.value()) if hasattr(self, 'iso_cores') else "1",
            },
            "Simulation": {
                'cores': str(self.sim_cores.value()) if hasattr(self, 'sim_cores') else "1",
                'batch_size': str(self.sim_batch_size.value()) if hasattr(self, 'sim_batch_size') else "1",
                'sample_perc': str(self.sim_sample_perc.value()) if hasattr(self, 'sim_sample_perc') else "10",
                'sample_type': self.sim_sample_type.currentText().lower() if hasattr(self,
                                                                                     'sim_sample_type') else "random",
                'use_isomorphism': self.use_isomorphism.isChecked() if hasattr(self, 'use_isomorphism') else "False",
                'energy_threshold': str(self.energy_threshold.value()) if hasattr(self, 'energy_threshold') else "20",
                'proton_numb': str(self.proton_numb.value()) if hasattr(self, 'proton_numb') else "1",
                'I': str(self.i_s.value()) if hasattr(self, 'i_s') else "0.25",
                'temp': str(self.temp.value()) if hasattr(self, 'temp') else "298.15",
                'min_pH': str(self.sim_min_ph.value()) if hasattr(self, 'sim_min_ph') else "0",
                'max_pH': str(self.sim_max_ph.value()) if hasattr(self, 'sim_max_ph') else "35",
                'step_pH': str(self.sim_step_ph.value()) if hasattr(self, 'sim_step_ph') else "0.5",
            },
            "CRN": {
                'Full_CRN': self.full_crn.currentText() if hasattr(self, 'full_crn') else "False",
                'Selected_model': str(self.selected_model.value()) if hasattr(self, 'selected_model') else "0",
                'Plot_3D': self.plot_3d.currentText() if hasattr(self, 'plot_3d') else "False",
            },
            "Scaling": {
                'scaling_mode': self.scaling_mode.currentText() if hasattr(self, 'scaling_mode') else "best_rmse",
                'experimental_set': self.exp_set.currentText() if hasattr(self, 'exp_set') else "",
            },
            "Speciation": {
                'labels_file': self.labels_file.text() if hasattr(self, 'labels_file') else "",
                'speciation_labels': self.spec_labels.text() if hasattr(self, 'spec_labels') else "all",
                'min_pH': str(self.spec_min_ph.value()) if hasattr(self, 'spec_min_ph') else "0",
                'max_pH': str(self.spec_max_ph.value()) if hasattr(self, 'spec_max_ph') else "14",
                'step_pH': str(self.spec_step_ph.value()) if hasattr(self, 'spec_step_ph') else "0.1",
                'cores': str(self.spec_cores.value()) if hasattr(self, 'spec_cores') else "1",
                'batch_size': str(self.spec_batch_size.value()) if hasattr(self, 'spec_batch_size') else "1",
                'm_idx': self.m_idx.currentText() if hasattr(self, 'm_idx') else "0",
                'npz_file': self.npz_file.text() if hasattr(self, 'npz_file') else "",
                'phase_dir': self.phase_dir.text() if hasattr(self, 'phase_dir') else "phase_diagram_%s" % system_name,
                'model_subset_file': self.model_subset_file.text() if hasattr(self, 'model_subset_file') else "",
            },
            "InternalConditions": {
                'restrain_addition': str(self.rest_add.value()) if hasattr(self, 'rest_add') else "1",
                'restrain_condensation': str(self.rest_cond.value()) if hasattr(self, 'rest_cond') else "0.0",
                'include_dimerization': self.include_dimerization.isChecked() if hasattr(self,
                                                                                         'include_dimerization') else "False",
                'force_stoich': self.force_sto.text() if hasattr(self, 'force_sto') else "0",
                'adjust_protons_hydration': self.adj_prot.isChecked() if hasattr(self, 'adj_prot') else "False",
            },
            "Clustering": {
                'cluster_dir': self.cluster_dir.text() if hasattr(self, 'cluster_dir') else "",
                'npz_cluster_file': self.npz_cluster_file.text() if hasattr(self, 'npz_cluster_file') else "",
                'features_file': self.features_file.text() if hasattr(self, 'features_file') else "",
                'n_clusters': str(self.n_clusters.value()) if hasattr(self, 'n_clusters') else "0",
                'normalize_feats': self.normalize_feats.currentText() if hasattr(self, 'normalize_feats') else "False",
            },
            "Visualization": {
                'col_dict': (
                    "" if (hasattr(self, 'clustering_col_dict') and
                           self.clustering_col_dict.col_dict_combo.currentText() == "None (Default)")
                    else self.clustering_col_dict.col_dict_combo.currentText() if hasattr(self, 'clustering_col_dict')
                    else "None"
                ),
                'plot_list': self.plot_list.text() if hasattr(self, 'plot_list') else "all",
                'boxplot_list': self.boxplot_list.text() if hasattr(self, 'boxplot_list') else "all",
            },
        }

        # Reference types (checkbox dict)
        if hasattr(self, 'ref_types') and isinstance(self.ref_types, dict):
            checked_types = [key for key, checkbox in self.ref_types.items() if checkbox.isChecked()]
            config_dict["Simulation"]["reference_types"] = ",".join(checked_types)

        # Feature selection (checkbox dict)
        if hasattr(self, 'sel_features') and isinstance(self.sel_features, dict):
            checked_types = [key for key, checkbox in self.sel_features.items() if checkbox.isChecked()]
            config_dict["Clustering"]["feats_list"] = ",".join(checked_types)

        # Cluster group selection (checkbox dict)
        if hasattr(self, 'selected_group_idx') and isinstance(self.selected_group_idx, dict):
            checked_groups = [key for key, checkbox in self.selected_group_idx.items() if checkbox.isChecked()]
            config_dict["Clustering"]["sel_groups"] = ",".join(checked_groups)

        # Mode-specific concentration and reference compound fields
        if is_ipa:
            config_dict['Simulation'].update({
                'C0': str(self.sim_c0.value()) if hasattr(self, 'sim_c0') else "0.1",
                'ref_compound': self.ref_compound.text() if hasattr(self, 'ref_compound') else "",
            })
            config_dict['Speciation'].update({
                'C': str(self.spec_C0.value()) if hasattr(self, 'spec_C0') else "0.1",
                'min_logC': str(self.phase_min_c.value()) if hasattr(self, 'phase_min_c') else "0",
                'max_logC': str(self.phase_max_c.value()) if hasattr(self, 'phase_max_c') else "14",
                'num_logC': str(self.phase_num_c.value()) if hasattr(self, 'phase_num_c') else "0.1",
            })
        else:
            config_dict['Simulation'].update({
                'CM': str(self.sim_CM.value()) if hasattr(self, 'sim_CM') else "0.1",
                'CX': str(self.sim_CX.value()) if hasattr(self, 'sim_CX') else "0.1",
                'ref_compound': (
                    f"{self.sim_M_ref_compound.text()},{self.sim_X_ref_compound.text()}"
                    if hasattr(self, 'sim_M_ref_compound') and hasattr(self, 'sim_X_ref_compound')
                    else ""
                ),
            })
            config_dict['Speciation'].update({
                'C_M': str(self.spec_CM.value()) if hasattr(self, 'spec_CM') else "0.1",
                'C_X': str(self.spec_CX.value()) if hasattr(self, 'spec_CX') else "0.1",
                'min_Ratio': str(self.min_metal_ratio.value()) if hasattr(self, 'min_metal_ratio') else "0.0001",
                'max_Ratio': str(self.max_metal_ratio.value()) if hasattr(self, 'max_metal_ratio') else "10.0",
                'num_Ratio': str(self.num_metal_ratio.value()) if hasattr(self, 'num_metal_ratio') else "2",
            })

        return config_dict

    def get_gui_params(self):
        """
        Return the current GUI parameters as a nested dictionary.

        Delegates to ``collect_config()`` which is the single source of truth
        for reading the GUI state.
        """
        return self.collect_config()

    def transform_absolute_to_relative_path(self, absolute_path):
        """
        Transform absolute paths to relative paths based on the project structure.
        Leave relative paths unchanged.

        Args:
            absolute_path (str): The path value to transform

        Returns:
            str: The transformed path (relative if it was absolute, unchanged if already relative)
        """
        if not absolute_path or not isinstance(absolute_path, str):
            return absolute_path

        # Check if path is already relative (doesn't start with / or drive letter on Windows)
        if not os.path.isabs(absolute_path):
            return absolute_path

        try:
            # Get the project root directory (parent of pomsimulator)
            project_root = os.path.dirname(script_dir) + "/pomsimulator"
            # Check if the absolute path is within inputs directory
            inputs_dir = os.path.join(project_root, "inputs")
            if absolute_path.startswith(inputs_dir):
                rel_path = os.path.relpath(absolute_path, inputs_dir)
                return "../inputs/" + rel_path.replace("\\", "/")

            # Check if the absolute path is within outputs directory
            outputs_dir = os.path.join(project_root, "outputs")
            if absolute_path.startswith(outputs_dir):
                rel_path = os.path.relpath(absolute_path, outputs_dir)
                return "../outputs/" + rel_path.replace("\\", "/")

            # For other paths within the project, make relative to project root
            try:
                rel_path = os.path.relpath(absolute_path, project_root)
                if not rel_path.startswith(".."):
                    return "../" + rel_path.replace("\\", "/")
                else:
                    # Path is outside project, keep as absolute
                    return absolute_path
            except ValueError:
                # Paths are on different drives (Windows), keep absolute
                return absolute_path

        except Exception as e:
            print(f"Error transforming path {absolute_path}: {e}")
            return absolute_path

    def browse_directory(self, line_edit):
        """
        Open a directory selection dialog and set the selected directory path to the provided line edit.
        Automatically transforms absolute paths to relative paths based on project structure.

        This method displays a directory dialog that allows the user to navigate and select
        a directory. The dialog starts in the pomsimulator directory within the script's
        directory. If the user selects a directory, its path is transformed to relative format
        and set as the text of the provided line edit widget.

        Parameters:
            line_edit (QLineEdit): The line edit widget where the selected directory path
                          will be displayed. This widget's text will be updated
                          with the relative path of the selected directory.

        Returns:
            None: This method updates the line_edit widget in-place and doesn't return a value.
        """
        # Define the pomsimulator root directory
        pomsimulator_root = os.path.join(script_dir)

        # Make sure the directory exists
        if not os.path.exists(pomsimulator_root):
            pomsimulator_root = script_dir

        # Open directory dialog starting from pomsimulator root
        directory = QFileDialog.getExistingDirectory(
            self,
            "Select Directory",
            pomsimulator_root,
            QFileDialog.DontUseNativeDialog
        )

        if directory:
            # Transform absolute path to relative path
            relative_path = self.transform_absolute_to_relative_path(directory)

            # Find which data key this line_edit is registered with
            for key, elements in self.ui_elements.items():
                if line_edit in elements:
                    # Update the shared data, which will update all related UI elements
                    self.update_shared_data(key, relative_path)
                    return

            # If not found in registered elements, just update the line edit directly
            line_edit.setText(relative_path)

    def browse_file(self, line_edit, filter_str="All Files (*)"):
        """
        Open a file selection dialog and set the selected file path to the provided line edit.
        Automatically transforms absolute paths to relative paths based on project structure.

        This method displays a file dialog that allows the user to navigate and select
        a file. The dialog starts in the pomsimulator directory within the script's
        directory. If the user selects a file, its path is transformed to relative format
        and set as the text of the provided line edit widget.

        Parameters:
            line_edit (QLineEdit): The line edit widget where the selected file path
                          will be displayed. This widget's text will be updated
                          with the relative path of the selected file.
            filter_str (str): Optional file filter string to limit the types of files
                     shown in the dialog. Defaults to "All Files (*)".
                     Example formats: "Text Files (*.txt)" or
                     "Images (*.png *.jpg)".

        Returns:
            None: This method updates the line_edit widget in-place and doesn't return a value.
        """
        # Define the pomsimulator root directory
        pomsimulator_root = os.path.join(script_dir)

        # Make sure the directory exists
        if not os.path.exists(pomsimulator_root):
            pomsimulator_root = script_dir

        # Open file dialog starting from pomsimulator root
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select File",
            pomsimulator_root,
            filter_str,
            options=QFileDialog.DontUseNativeDialog
        )

        if file_path:
            # Transform absolute path to relative path
            relative_path = self.transform_absolute_to_relative_path(file_path)
            line_edit.setText(relative_path)

    def create_mode_toolbar(self):
        """
        Create a toolbar with IPA/HPA mode selection buttons.
        """
        try:
            # Create toolbar
            self.mode_toolbar = self.addToolBar("Mode Selection")
            self.mode_toolbar.setMovable(False)  # Prevent moving the toolbar

            # Create button group for exclusive selection
            self.mode_button_group = QButtonGroup()

            # Create IPA button
            self.ipa_mode_btn = QPushButton("IPA")
            self.ipa_mode_btn.setCheckable(True)
            self.ipa_mode_btn.setToolTip("Switch to Isopolyanion mode")
            self.ipa_mode_btn.setStyleSheet("""
                QPushButton {
                    font-size: 12pt;
                    font-weight: bold;
                    padding: 8px 16px;
                    margin: 2px;
                    border: 2px solid #cccccc;
                    border-radius: 5px;
                    background-color: #f4f4f4;
                    color: black;
                }
                QPushButton:checked {
                    background-color: #4CAF50;
                    color: black;
                    border-color: #45a049;
                }
                QPushButton:hover {
                    background-color: #e8e8e8;
                    color: black;
                }
                QPushButton:checked:hover {
                    background-color: #45a049;
                    color: black;
                }
            """)
            # Create HPA button
            self.hpa_mode_btn = QPushButton("HPA")
            self.hpa_mode_btn.setCheckable(True)
            self.hpa_mode_btn.setToolTip("Switch to Heteropolyanion mode")
            self.hpa_mode_btn.setStyleSheet("""
                QPushButton {
                    font-size: 12pt;
                    font-weight: bold;
                    padding: 8px 16px;
                    margin: 2px;
                    border: 2px solid #cccccc;
                    border-radius: 5px;
                    background-color: #f4f4f4;
                    color: black;
                }
                QPushButton:checked {
                    background-color: #2196F3;
                    color: black;
                    border-color: #1976D2;
                }
                QPushButton:hover {
                    background-color: #e8e8e8;
                    color: black;
                }
                QPushButton:checked:hover {
                    background-color: #1976D2;
                    color: black;
                }
            """)

            # Add buttons to button group for exclusive selection
            self.mode_button_group.addButton(self.ipa_mode_btn, 0)
            self.mode_button_group.addButton(self.hpa_mode_btn, 1)

            # Connect signals
            self.ipa_mode_btn.clicked.connect(lambda: self.request_mode_change("IPA"))
            self.hpa_mode_btn.clicked.connect(lambda: self.request_mode_change("HPA"))

            # Add buttons to toolbar
            self.mode_toolbar.addWidget(self.ipa_mode_btn)
            self.mode_toolbar.addWidget(self.hpa_mode_btn)

            # Add separator and mode label
            self.mode_toolbar.addSeparator()
            self.mode_label = QLabel("Current Mode: IPA")
            self.mode_label.setStyleSheet("font-weight: bold; margin: 5px;")
            self.mode_toolbar.addWidget(self.mode_label)

            # Set initial mode
            self.global_mode = "IPA"
            self.ipa_mode_btn.setChecked(True)

        except Exception as e:
            print(f"Error creating mode toolbar: {e}")

    def toggle_mode_with_shortcut(self):
        """
        Toggle mode with keyboard shortcut.
        """
        try:
            current_mode = getattr(self, 'global_mode', 'IPA')
            new_mode = "HPA" if current_mode == "IPA" else "IPA"

            self.request_mode_change(new_mode)
        except Exception as e:
            print(f"Error in toggle_mode_with_shortcut: {e}")

    def refresh_interface(self):
        """
        Refresh the interface with new mode.
        """
        try:
            reply = QMessageBox.question(
                self,
                "Confirm GUI refresh",
                "Are you sure you want to refresh the interface?\n\n"
                "This will reset all current parameters and reload the interface.",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )
            if reply == QMessageBox.Yes:
                self.update_console("Refreshing interface ...")
                self.clear_all_data_and_ui()
                self.refresh_all_ui_elements()
                self.repaint()
                self.clear_console()
                self.update_console("Interface refreshed.")
            else:
                # User cancelled - revert button selection
                self.update_console("Interface refresh cancelled.")

        except Exception as e:
            print(f"Error in refresh_interface: {e}")

    def request_mode_change(self, new_mode):
        """
        Request a mode change with user confirmation.

        Args:
            new_mode (str): The requested new mode ("IPA" or "HPA")
        """
        try:
            # If it's the same mode, do nothing
            if hasattr(self, 'global_mode') and self.global_mode == new_mode:
                return

            # Show confirmation dialog
            mode_name = "Isopolyoxometalates" if new_mode == "IPA" else "Heteropolyoxometalates"
            current_mode_name = "Isopolyoxometalates" if getattr(self, 'global_mode',
                                                                 'IPA') == "IPA" else "Heteropolyoxometalates"

            reply = QMessageBox.question(
                self,
                "Confirm Mode Change",
                f"Are you sure you want to switch from {current_mode_name} to {mode_name} mode?\n\n"
                f"This will reset all current parameters and reload the interface.",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No
            )

            if reply == QMessageBox.Yes:
                # User confirmed - proceed with mode change
                self.execute_mode_change(new_mode)
            else:
                # User cancelled - revert button selection
                self.revert_mode_button_selection()

        except Exception as e:
            print(f"Error in request_mode_change: {e}")
            self.revert_mode_button_selection()

    def revert_mode_button_selection(self):
        """
        Revert the mode button selection to the current global mode.
        """
        try:
            current_mode = getattr(self, 'global_mode', 'IPA')

            # Block signals to prevent triggering another mode change
            self.ipa_mode_btn.blockSignals(True)
            self.hpa_mode_btn.blockSignals(True)

            # Set the correct button as checked
            if current_mode == "IPA":
                self.ipa_mode_btn.setChecked(True)
                self.hpa_mode_btn.setChecked(False)
            elif current_mode == "HPA":
                self.hpa_mode_btn.setChecked(True)
                self.ipa_mode_btn.setChecked(False)

            # Unblock signals
            self.ipa_mode_btn.blockSignals(False)
            self.hpa_mode_btn.blockSignals(False)

        except Exception as e:
            print(f"Error reverting mode button selection: {e}")

    def execute_mode_change(self, new_mode):
        """
        Execute the actual mode change after user confirmation.

        Args:
            new_mode (str): The new mode to switch to ("IPA" or "HPA")
        """
        try:
            old_mode = getattr(self, 'global_mode', 'IPA')

            # Update global mode
            self.global_mode = new_mode
            self.pom_type = new_mode  # Keep compatibility with existing code

            # Update toolbar button selection
            if new_mode == "IPA":
                self.ipa_mode_btn.setChecked(True)
                self.hpa_mode_btn.setChecked(False)
            elif new_mode == "HPA":
                self.hpa_mode_btn.setChecked(True)
                self.ipa_mode_btn.setChecked(False)

            # Update mode label
            mode_name = "Isopolyoxometalates" if new_mode == "IPA" else "Heteropolyoxometalates"
            self.mode_label.setText(f"Current Mode: {new_mode}")

            # Clear console and add mode change message
            self.clear_console()
            self.update_console(f"Mode changed from {old_mode} to {new_mode} - rebuilding interface", "system")

            # Clear all data and UI values
            self.clear_all_data_and_ui()

            # Rebuild mode-dependent tabs
            self.rebuild_interface_for_mode(new_mode)

            # Reset cached data and button states
            self.reset_cached_data()
            self.reset_all_button_states()
            self.reset_all_tab_status()
            # Update window title
            self.update_window_title_for_mode(new_mode)

            self.repaint()  # Trigger a repaint to update the interface

            self.update_console(f"Interface rebuilt for {mode_name} mode", "system")

            QTimer.singleShot(2000, lambda: self.clear_console())

        except Exception as e:
            print(f"Error executing mode change: {e}")
            self.revert_mode_button_selection()

    def clear_all_data_and_ui(self):
        """
        Clear all shared data and UI element values.
        This is the main clearing function that handles both data and visual elements.
        """
        try:
            self.update_console("Clearing all data and UI values...", "system")

            # Clear shared data
            if hasattr(self, 'shared_data'):
                self.shared_data.clear()

            # Clear stored configuration data
            if hasattr(self, 'current_config'):
                self.current_config = None

            # Clear all UI element values visually
            self.clear_all_ui_element_values()

            self.update_console("All data and UI values cleared successfully", "system")

        except Exception as e:
            print(f"Error clearing all data and UI: {e}")
            self.update_console(f"Error clearing all data and UI: {e}", "system")

    def clear_all_ui_element_values(self):
        """
        Clear all registered UI element values to their defaults.
        This preserves the UI elements registry but resets their visual values.
        """
        try:
            # Iterate through all registered UI elements and clear their values
            for key, elements in self.ui_elements.items():
                for element in elements:
                    if self.is_widget_valid(element):
                        self.clear_ui_element_value(element)

        except Exception as e:
            print(f"Error clearing UI element values: {e}")

    def clear_ui_element_value(self, element):
        """
        Clear the value of a single UI element based on its type.

        Parameters:
            element (QWidget): The UI element to clear
        """
        try:
            if isinstance(element, QLineEdit):
                element.clear()
            elif isinstance(element, QSpinBox):
                element.setValue(0)
            elif isinstance(element, CustomDoubleSpinBox):
                element.setValue(0.0)
            elif isinstance(element, QComboBox):
                element.setCurrentIndex(0)
            elif isinstance(element, QCheckBox):
                element.setChecked(False)
            elif isinstance(element, QTextEdit):
                element.clear()
            # Add more widget types as needed
        except Exception as e:
            print(f"Error clearing value for element {element}: {e}")

    def rebuild_interface_for_mode(self, mode):
        """
        Rebuild the entire interface for the specified mode.

        Args:
            mode (str): The mode to rebuild for ("IPA" or "HPA")
        """
        try:
            # Update all tabs for the new mode
            self.update_all_tabs_for_mode(mode)

            # Force a complete repaint
            self.repaint()

        except Exception as e:
            print(f"Error rebuilding interface for mode {mode}: {e}")

    def update_all_tabs_for_mode(self, mode):
        """
        Update all tabs for the new mode.

        Args:
            mode (str): The new mode ("IPA" or "HPA")
        """
        try:
            # Update mode-dependent tabs (these need to be rebuilt for the new mode)
            if hasattr(self, 'sim_params_layout'):
                self.update_sim_parameters()
            if hasattr(self, 'speciation_params_layout'):
                self.update_speciation_parameters()
            if hasattr(self, 'phase_params_layout'):
                self.update_phase_parameters()
            if hasattr(self, 'plot_spec_params_layout'):
                self.update_plot_speciation_parameters()
            if hasattr(self, 'plot_phase_params_layout'):
                self.update_plot_phase_parameters()

            # After rebuilding mode-dependent tabs, refresh all UI elements
            self.refresh_all_ui_elements()

        except Exception as e:
            print(f"Error updating tabs for mode {mode}: {e}")

    def refresh_all_ui_elements(self):
        """
        Refresh all UI elements by updating them from the shared data.
        This ensures all tabs reflect the current values.
        """
        try:
            if not hasattr(self, 'shared_data') or not hasattr(self, 'ui_elements'):
                return

            # Iterate through all registered UI elements and update them
            for data_key, elements in self.ui_elements.items():
                if data_key in self.shared_data:
                    value = self.shared_data[data_key]

                    # Update all elements for this data key
                    for element in elements:
                        if self.is_widget_valid(element):
                            try:
                                self.set_ui_element_value(element, value)
                            except Exception as e:
                                print(f"Error updating element for {data_key}: {e}")

            self.update_console("All UI elements refreshed", "system")

        except Exception as e:
            print(f"Error refreshing UI elements: {e}")

    def reset_cached_data(self):
        """
        Reset any cached data or temporary values.
        """
        try:
            # Reset any runners
            if hasattr(self, 'mol_file_runner'):
                self.mol_file_runner = None
            if hasattr(self, 'isomorph_runner'):
                self.isomorph_runner = None
            if hasattr(self, 'simulation_runner'):
                self.simulation_runner = None
            if hasattr(self, 'speciation_runner'):
                self.speciation_runner = None
            # TODO implement any other cached data specific to your application

            # Reset any other cached data specific to your application
            # Add more resets as needed based on your application's state

        except Exception as e:
            print(f"Error resetting cached data: {e}")

    def reset_all_button_states(self):
        """
        Reset all button states to their default values.
        """
        try:
            # Reset presimulation buttons
            if hasattr(self, 'gen_mol_run_btn'):
                self.gen_mol_run_btn.setEnabled(True)
            if hasattr(self, 'gen_mol_stop_btn'):
                self.gen_mol_stop_btn.setEnabled(False)
            if hasattr(self, 'comp_iso_run_btn'):
                self.comp_iso_run_btn.setEnabled(True)
            if hasattr(self, 'comp_iso_stop_btn'):
                self.comp_iso_stop_btn.setEnabled(False)

            # Reset simulation buttons
            if hasattr(self, 'ipa_sim_run_btn'):
                self.ipa_sim_run_btn.setEnabled(True)
            if hasattr(self, 'ipa_sim_stop_btn'):
                self.ipa_sim_stop_btn.setEnabled(False)
            if hasattr(self, 'hpa_sim_run_btn'):
                self.hpa_sim_run_btn.setEnabled(True)
            if hasattr(self, 'hpa_sim_stop_btn'):
                self.hpa_sim_stop_btn.setEnabled(False)

            # Reset speciation buttons
            if hasattr(self, 'ipa_comp_spec_run_btn'):
                self.ipa_comp_spec_run_btn.setEnabled(True)
            if hasattr(self, 'ipa_comp_spec_stop_btn'):
                self.ipa_comp_spec_stop_btn.setEnabled(False)
            if hasattr(self, 'hpa_comp_spec_run_btn'):
                self.hpa_comp_spec_run_btn.setEnabled(True)
            if hasattr(self, 'hpa_comp_spec_stop_btn'):
                self.hpa_comp_spec_stop_btn.setEnabled(False)

        except Exception as e:
            print(f"Error resetting button states: {e}")

    def reset_all_tab_status(self):
        for tab in self.tab_status:
            if type(self.tab_status[tab]) == dict:
                for subtab in self.tab_status[tab]:
                    self.set_tab_idle(tab, subtab)
            else:
                subtab = None
                self.set_tab_idle(tab, subtab)

    def update_window_title_for_mode(self, mode):
        """
        Update the window title to reflect the current mode.

        Args:
            mode (str): The current mode ("IPA" or "HPA")
        """
        try:
            base_title = "POMSimulator"
            mode_name = "Isopolyoxometalates" if mode == "IPA" else "Heteropolyoxometalates"
            new_title = f"{base_title} - {mode_name} Mode"
            self.setWindowTitle(new_title)
        except Exception as e:
            print(f"Error updating window title: {e}")

    def is_widget_valid(self, widget):
        """
        Check if a widget is still valid and not deleted.

        Args:
            widget: The widget to check

        Returns:
            bool: True if the widget is valid, False otherwise
        """
        try:
            # Try to access a basic property of the widget
            if widget is None:
                return False

            # For Qt widgets, try to access the parent property
            # This will raise RuntimeError if the widget has been deleted
            _ = widget.parent()

            # Additional check: make sure it's actually a widget
            if not hasattr(widget, 'isVisible'):
                return False

            return True

        except (RuntimeError, AttributeError):
            # Widget has been deleted or is not a valid Qt widget
            return False
        except Exception as e:
            print(f"Unexpected error checking widget validity: {e}")
            return False

    def register_ui_element(self, element, data_key):
        """
        Register a UI element to be synchronized with shared data.

        Parameters:
            element: The UI element to register
            data_key: The key in the shared_data dictionary
        """
        # Create the list for this key if it doesn't exist
        if data_key not in self.ui_elements:
            self.ui_elements[data_key] = []
        elif not isinstance(self.ui_elements[data_key], list):
            # Convert single widget to list (from get_or_create_widget)
            existing_widget = self.ui_elements[data_key]
            self.ui_elements[data_key] = [existing_widget]

        # Add the element to the list if not already present
        if element not in self.ui_elements[data_key]:
            self.ui_elements[data_key].append(element)

        # Disconnect any existing connections to avoid duplicates
        try:
            if isinstance(element, QLineEdit):
                try:
                    element.textChanged.disconnect()
                except TypeError:
                    pass
                element.textChanged.connect(lambda text, e=element, k=data_key: self.sync_ui_elements(k, text, e))
            elif isinstance(element, QSpinBox):
                try:
                    element.valueChanged.disconnect()
                except TypeError:
                    pass
                element.valueChanged.connect(lambda value, e=element, k=data_key: self.sync_ui_elements(k, value, e))
            elif isinstance(element, QDoubleSpinBox):
                try:
                    element.valueChanged.disconnect()
                except TypeError:
                    pass
                element.valueChanged.connect(lambda value, e=element, k=data_key: self.sync_ui_elements(k, value, e))
            elif isinstance(element, QComboBox):
                try:
                    element.currentTextChanged.disconnect()
                except TypeError:
                    pass
                element.currentTextChanged.connect(
                    lambda text, e=element, k=data_key: self.sync_ui_elements(k, text, e))
            elif isinstance(element, QCheckBox):
                try:
                    element.stateChanged.disconnect()
                except TypeError:
                    pass
                element.stateChanged.connect(
                    lambda state, e=element, k=data_key: self.sync_ui_elements(k, state == Qt.Checked, e))
        except Exception as e:
            print(f"Error connecting signal for {data_key}: {str(e)}")

    def sync_ui_elements(self, data_key, value, sender=None):
        """
        Synchronize all UI elements that share the same data key.

        Parameters:
            data_key: The key in the shared_data dictionary
            value: The new value to set
            sender: The UI element that triggered the update (to avoid loops)
        """
        # Update the shared data
        self.shared_data[data_key] = value

        # Update all UI elements with this data key, except the sender
        if data_key in self.ui_elements:
            # Create a copy of the list to avoid modification during iteration
            elements_copy = list(self.ui_elements[data_key])

            for element in elements_copy:
                if element != sender:  # Skip the sender to avoid loops
                    try:
                        # Check if the widget is still valid
                        if not self.is_widget_valid(element):
                            # Remove invalid widget from the list
                            if element in self.ui_elements[data_key]:
                                self.ui_elements[data_key].remove(element)
                            continue

                        if isinstance(element, QLineEdit):
                            element.setText(str(value))
                        elif isinstance(element, CustomDoubleSpinBox):
                            element.setValue(float(value))
                        elif isinstance(element, QSpinBox):
                            element.setValue(int(value))
                        elif isinstance(element, QComboBox):
                            index = element.findText(str(value))
                            if index >= 0:
                                element.setCurrentIndex(index)
                        elif isinstance(element, QCheckBox):
                            element.setChecked(bool(value))
                    except RuntimeError as e:
                        # Handle "wrapped C/C++ object has been deleted" error
                        print(f"Widget deleted for {data_key}, removing from list: {e}")
                        if element in self.ui_elements[data_key]:
                            self.ui_elements[data_key].remove(element)
                    except Exception as e:
                        print(f"Error updating {data_key} for element {element}: {str(e)}")
                        # Remove problematic element from the list
                        if element in self.ui_elements[data_key]:
                            self.ui_elements[data_key].remove(element)

    def update_shared_data(self, key, value, source_element=None):
        """
        Update shared data and synchronize all UI elements.

        Parameters:
            key (str): The key in the shared_data dictionary to update
            value: The new value to set
            source_element (QWidget): The UI element that triggered the update (to avoid loops)
        """
        # Update the shared data
        self.shared_data[key] = value

        # Update all registered UI elements except the source
        for element in self.ui_elements.get(key, []):
            if element != source_element:
                self.set_ui_element_value(element, value)

    def set_ui_element_value(self, element, value):
        """
        Set the value of a UI element based on its type.

        Parameters:
            element (QWidget): The UI element to update
            value: The value to set
        """
        try:
            if isinstance(element, QLineEdit):
                element.setText(str(value))
            elif isinstance(element, QSpinBox):
                element.setValue(int(value))
            elif isinstance(element, CustomDoubleSpinBox):
                element.setValue(float(value))
            elif isinstance(element, QComboBox):
                index = element.findText(str(value))
                if index >= 0:
                    element.setCurrentIndex(index)
            # Add more widget types as needed
        except (ValueError, TypeError) as e:
            print(f"Error setting value {value} to element {element}: {e}")

    def update_console(self, message, tab_name=None):
        """
        Update the shared console with a new message.

        Parameters:
            message (str): The message to add to the console
            tab_name (str, optional): The name of the tab (used for context in the message).
        """
        try:
            # Add timestamp to message
            from datetime import datetime
            timestamp = datetime.now().strftime("%H:%M:%S")

            # Add tab prefix if tab_name is provided for context
            if tab_name:
                formatted_message = f"[{timestamp}] [{tab_name.upper()}] {message}"
            else:
                formatted_message = f"[{timestamp}] {message}"

            # Check if shared console exists
            if hasattr(self, 'shared_console') and self.shared_console is not None:
                # Make sure it's actually a QTextEdit widget
                if hasattr(self.shared_console, 'append'):
                    self.shared_console.append(formatted_message)
                    # Auto-scroll to bottom
                    if hasattr(self.shared_console, 'verticalScrollBar'):
                        scrollbar = self.shared_console.verticalScrollBar()
                        scrollbar.setValue(scrollbar.maximum())
                else:
                    print(f"Shared console object doesn't have append method. Type: {type(self.shared_console)}")
            else:
                print(f"Shared console not found. Fallback: {message}")
                # Fallback to print if console fails
                print(f"CONSOLE: {message}")

        except Exception as e:
            print(f"Console error: {e}, Message: {message}")
            # Fallback to print if console fails
            print(f"CONSOLE: {message}")

    def add_shared_console(self, layout):
        """
        Add the shared console to the main window layout.

        This creates a single console that all tabs will use for displaying output.
        The console is placed at the bottom of the main window and is always visible.

        Parameters:
            layout (QVBoxLayout): The main window layout where the console will be added.
        """
        console_group = QGroupBox("Console Output")
        console_layout = QVBoxLayout()

        # Create console controls
        controls_layout = QHBoxLayout()
        clear_btn = QPushButton("Clear")
        clear_btn.setMaximumWidth(80)
        clear_btn.clicked.connect(self.clear_console)
        controls_layout.addWidget(clear_btn)
        # controls_layout.addStretch()  # Push clear button to the left

        console_layout.addLayout(controls_layout)

        # Create the shared console
        self.shared_console = QTextEdit()
        self.shared_console.setReadOnly(True)
        # self.shared_console.setMaximumHeight(200)  # Limit height so it doesn't dominate the UI
        self.shared_console.setMinimumHeight(200)  # Ensure it's always visible

        # Set monospace font for better readability
        monospace_font = QFont("Monospace")
        monospace_font.setStyleHint(QFont.Monospace)
        monospace_font.setPointSize(11)
        self.shared_console.setFont(monospace_font)

        # Apply initial theme-based styling
        self.update_console_theme()

        console_layout.addWidget(self.shared_console)
        console_group.setLayout(console_layout)
        layout.addWidget(console_group)

    def update_console_theme(self):
        """
        Apply the current theme-based styling to the shared console.
        """
        if not hasattr(self, 'shared_console') or self.shared_console is None:
            return
        current_theme = getattr(self, 'current_theme', 'light')

        if self.dark_theme:
            self.shared_console.setStyleSheet("""
                QTextEdit {
                    background-color: #1e1e1e;
                    color: #ffffff;
                    border: 1px solid #555555;
                    border-radius: 4px;
                    selection-background-color: #3d3d3d;
                    selection-color: #ffffff;
                }
            """)
        else:
            self.shared_console.setStyleSheet("""
                        QTextEdit {
                            background-color: #f8f8f8;
                            color: #000000;
                            border: 1px solid #cccccc;
                            border-radius: 4px;
                            selection-background-color: #316AC5;
                            selection-color: #ffffff;
                        }
                    """)

    def clear_console(self):
        """
        Clear the shared console output.

        Parameters:
            tab_name (str, optional): Kept for compatibility, but ignored since
                                     we now use a single shared console.
        """
        if hasattr(self, 'shared_console') and self.shared_console is not None:
            self.shared_console.clear()
            self.update_console("Console cleared")
        else:
            print("No shared console to clear")

    """Workers"""

    def _get_runner_registry(self):
        """Return a list of runner descriptors for all known runners.

        Each descriptor is a dict with:
            runner_attr  (str)  – attribute name on self that holds the runner
            run_btn_attr (str or callable) – attribute name of the run button, or a
                                             zero-argument callable that returns it
            stop_btn_attr (str or callable) – attribute name of the stop button, or a
                                              zero-argument callable that returns it
            tab_name     (str)  – main tab name for set_tab_* calls
            subtab_name  (str or None) – subtab name for set_tab_* calls
            console_ctx  (str)  – context label for update_console calls
            stop_msg     (str)  – message to print when stopping

        For mode-dependent runners (simulation, speciation, phase), the button
        attributes are resolved lazily via callables so they always reflect the
        current pom_type.
        """

        def _sim_run():
            return self.ipa_comp_lgkf_run_btn if self.pom_type == "IPA" else self.hpa_comp_lgkf_run_btn

        def _sim_stop():
            return self.ipa_comp_lgkf_stop_btn if self.pom_type == "IPA" else self.hpa_comp_lgkf_stop_btn

        def _spec_run():
            return self.ipa_comp_spec_run_btn if self.pom_type == "IPA" else self.hpa_comp_spec_run_btn

        def _spec_stop():
            return self.ipa_comp_spec_stop_btn if self.pom_type == "IPA" else self.hpa_comp_spec_stop_btn

        def _phase_run():
            return self.ipa_phase_run_btn if self.pom_type == "IPA" else self.hpa_phase_run_btn

        def _phase_stop():
            return self.ipa_phase_stop_btn if self.pom_type == "IPA" else self.hpa_phase_stop_btn

        return [
            {
                'runner_attr': 'mol_file_runner',
                'run_btn_attr': 'gen_mol_run_btn',
                'stop_btn_attr': 'gen_mol_stop_btn',
                'tab_name': 'presimulation',
                'subtab_name': 'generate_mol',
                'console_ctx': 'presimulation',
                'stop_msg': 'Stopping mol file generation...',
            },
            {
                'runner_attr': 'isomorph_runner',
                'run_btn_attr': 'comp_iso_run_btn',
                'stop_btn_attr': 'comp_iso_stop_btn',
                'tab_name': 'presimulation',
                'subtab_name': 'compute_iso',
                'console_ctx': 'presimulation',
                'stop_msg': 'Stopping compute isomorphism...',
            },
            {
                'runner_attr': 'simulation_runner',
                'run_btn_attr': _sim_run,
                'stop_btn_attr': _sim_stop,
                'tab_name': 'simulation',
                'subtab_name': 'compute_lgkf',
                'console_ctx': 'simulation',
                'stop_msg': 'Stopping simulation...',
            },
            {
                'runner_attr': 'scaling_runner',
                'run_btn_attr': 'scale_run_btn',
                'stop_btn_attr': 'scale_stop_btn',
                'tab_name': 'scaling',
                'subtab_name': None,
                'console_ctx': 'scale_constants',
                'stop_msg': 'Stopping scaling...',
            },
            {
                'runner_attr': 'speciation_runner',
                'run_btn_attr': _spec_run,
                'stop_btn_attr': _spec_stop,
                'tab_name': 'speciation',
                'subtab_name': 'spec_diag',
                'console_ctx': 'speciation',
                'stop_msg': 'Stopping speciation...',
            },
            {
                'runner_attr': 'phase_runner',
                'run_btn_attr': _phase_run,
                'stop_btn_attr': _phase_stop,
                'tab_name': 'speciation',
                'subtab_name': 'phase_diag',
                'console_ctx': 'speciation',
                'stop_msg': 'Stopping phase identification...',
            },
            {
                'runner_attr': 'cluster_runner',
                'run_btn_attr': 'cluster_run_btn',
                'stop_btn_attr': 'cluster_stop_btn',
                'tab_name': 'clustering',
                'subtab_name': 'clustering',
                'console_ctx': 'clustering',
                'stop_msg': 'Stopping cluster analysis...',
            },
            {
                'runner_attr': 'selection_runner',
                'run_btn_attr': 'selection_run_btn',
                'stop_btn_attr': 'selection_stop_btn',
                'tab_name': 'clustering',
                'subtab_name': 'selection',
                'console_ctx': 'clustering',
                'stop_msg': 'Stopping selection...',
            },
            {
                'runner_attr': 'filtering_runner',
                'run_btn_attr': 'filtering_run_btn',
                'stop_btn_attr': 'filtering_stop_btn',
                'tab_name': 'clustering',
                'subtab_name': 'filtering',
                'console_ctx': 'clustering',
                'stop_msg': 'Stopping filtering...',
            },
            {
                'runner_attr': 'plot_spec_runner',
                'run_btn_attr': 'plot_spec_run_btn',
                'stop_btn_attr': 'plot_spec_stop_btn',
                'tab_name': 'plotting',
                'subtab_name': 'plot_spec',
                'console_ctx': 'plotting',
                'stop_msg': 'Stopping speciation plot...',
            },
            {
                'runner_attr': 'plot_phase_runner',
                'run_btn_attr': 'plot_phase_run_btn',
                'stop_btn_attr': 'plot_phase_stop_btn',
                'tab_name': 'plotting',
                'subtab_name': 'plot_phase',
                'console_ctx': 'plotting',
                'stop_msg': 'Stopping phase plot...',
            },
        ]

    def _resolve_btn(self, btn_attr):
        """Resolve a button attribute: either call it (if callable) or look it up on self."""
        if callable(btn_attr):
            try:
                return btn_attr()
            except AttributeError:
                return None
        return getattr(self, btn_attr, None)

    def _handle_runner_event(self, event, error_message=None):
        """Central dispatcher for runner lifecycle events.

        Parameters:
            event (str): One of 'finished', 'stopped', 'error'.
            error_message (str, optional): Error message for 'error' events.
        """
        sender = self.sender()
        for entry in self._get_runner_registry():
            runner = getattr(self, entry['runner_attr'], None)
            if runner is None or sender is not runner:
                continue

            run_btn = self._resolve_btn(entry['run_btn_attr'])
            stop_btn = self._resolve_btn(entry['stop_btn_attr'])

            if stop_btn is not None:
                stop_btn.setEnabled(False)
            if run_btn is not None:
                run_btn.setEnabled(True)

            tab_name = entry['tab_name']
            subtab_name = entry['subtab_name']
            ctx = entry['console_ctx']

            if event == 'finished':
                self.set_tab_success(tab_name, subtab_name)
                self.update_console("Function execution completed", ctx)
            elif event == 'stopped':
                self.set_tab_idle(tab_name, subtab_name)
                self.update_console("Function execution stopped by user", ctx)
            elif event == 'error':
                self.set_tab_error(tab_name, subtab_name)
                if error_message:
                    self.update_console(f"Function error: {error_message}", ctx)
            break  # Each sender matches at most one registry entry

    def stop_function(self):
        """Stop the currently running function."""
        for entry in self._get_runner_registry():
            runner = getattr(self, entry['runner_attr'], None)
            if runner is not None and runner.is_running:
                runner.stop()
                self.update_console(entry['stop_msg'], entry['console_ctx'])
                return

    def on_function_finished(self):
        """Handle function completion."""
        self._handle_runner_event('finished')

    def on_function_stopped(self):
        """Handle function stopped by user."""
        self._handle_runner_event('stopped')

    def on_function_error(self, error_message):
        """Handle function errors."""
        self._handle_runner_event('error', error_message=error_message)
        QMessageBox.critical(self, "Function Error", error_message)
        print(f"Function error: {error_message}")

    """Tabs"""

    def update_tab_status(self, tab_name, status, subtab_name=None):
        """
        Update the status indicator for a specific tab or subtab.

        Parameters:
            tab_name (str): Name of the main tab ('presimulation', 'simulation', etc.)
            status (str): Status to set ('idle', 'running', 'success', 'error')
            subtab_name (str, optional): Name of the subtab if updating a subtab
        """
        if tab_name not in self.tab_status:
            print(f"Tab '{tab_name}' not found in tab_status")
            return

        # Check if this tab has subtabs
        if isinstance(self.tab_status[tab_name], dict):
            # Tab has subtabs
            if subtab_name:
                if subtab_name in self.tab_status[tab_name]:
                    self.tab_status[tab_name][subtab_name] = status
            else:
                # No subtab specified, update all subtabs
                for subtab in self.tab_status[tab_name]:
                    self.tab_status[tab_name][subtab] = status
        else:
            # Tab has no subtabs
            if not subtab_name:
                self.tab_status[tab_name] = status

        # Refresh all indicators
        self.refresh_tab_indicators()

    def get_subtab_mapping(self, tab_key):
        """
        Get the mapping of subtab keys to their display names and indices.

        Parameters:
            tab_key (str): The main tab key ('presimulation', 'simulation', etc.)

        Returns:
            dict: Mapping of subtab_key -> (display_name, index)
        """
        subtab_mappings = {
            'presimulation': {
                'generate_mol': ('Generate Molfiles', 0),
                'compute_iso': ('Compute Isomorphism', 1)
            },
            'simulation': {
                'compute_lgkf': ('Compute lgkf', 0)
                # 'generate_crn': ('Generate CRN', 1)
            },
            'speciation': {
                'spec_diag': ('Speciation Diagram', 0),
                'phase_diag': ('Phase Diagram', 1)
            },
            'clustering': {
                'clustering': ('Clustering', 0),
                'selection': ('Model Selection', 1),
                'filtering': ('Boxplot Filtering', 2)
            },
            'plotting': {
                'plot_spec': ('Plot Speciation Diagram', 0),
                'plot_phase': ('Plot Phase Diagram', 1)
            }
        }

        return subtab_mappings.get(tab_key, {})

    def get_aggregate_tab_status(self, tab_key):
        """
        Get the aggregate status for a tab that contains subtabs.

        Priority order: error > running > success > idle
        """
        if tab_key not in self.tab_status:
            return 'idle'

        if isinstance(self.tab_status[tab_key], dict):
            # Tab has subtabs - calculate aggregate status
            subtab_statuses = list(self.tab_status[tab_key].values())

            # Priority order: error > running > success > idle
            if 'error' in subtab_statuses:
                return 'error'
            elif 'running' in subtab_statuses:
                return 'running'
            elif all(status == 'success' for status in subtab_statuses):
                return 'success'
            elif any(status == 'idle' for status in subtab_statuses):
                return 'idle'  # Partial success still shows as success
            else:
                return 'idle'
        else:
            # Tab has no subtabs - return direct status
            return self.tab_status[tab_key]

    def refresh_tab_indicators(self):
        """
        Refresh all tab status indicators based on current statuses.
        """
        status_indicators = {
            'idle': '∘',
            'running': '[♞]',
            'success': '[✓]',
            'error': '[✗]'
        }

        # Map tab names to their display names and indices
        tab_mapping = {
            'presimulation': ('Presimulation', 0),
            'simulation': ('Simulation', 1),
            'scaling': ('Scaling', 2),
            'speciation': ('Speciation', 3),
            'clustering': ('Clustering', 4),
            'plotting': ('Plotting', 5)
        }

        for tab_key, (display_name, index) in tab_mapping.items():
            if hasattr(self, 'tabs') and index < self.tabs.count():

                # Check if this tab has subtabs
                if isinstance(self.tab_status[tab_key], dict):
                    # Tab has subtabs - handle both subtabs and main tab

                    # First, update subtab indicators
                    main_tab_widget = self.tabs.widget(index)
                    if hasattr(main_tab_widget, 'findChild'):
                        subtab_widget = main_tab_widget.findChild(QTabWidget)
                        if subtab_widget:
                            # Update each subtab indicator
                            subtab_mapping = self.get_subtab_mapping(tab_key)
                            for subtab_key, (subtab_display_name, subtab_index) in subtab_mapping.items():
                                if subtab_index < subtab_widget.count():
                                    subtab_status = self.tab_status[tab_key].get(subtab_key, 'idle')
                                    subtab_indicator = status_indicators.get(subtab_status, '')

                                    if subtab_indicator:
                                        subtab_text = f"{subtab_display_name} {subtab_indicator}"
                                    else:
                                        subtab_text = subtab_display_name

                                    subtab_widget.setTabText(subtab_index, subtab_text)

                    # Then, calculate and set aggregate status for main tab
                    aggregate_status = self.get_aggregate_tab_status(tab_key)
                    main_indicator = status_indicators.get(aggregate_status, '')

                    if main_indicator:
                        main_tab_text = f"{display_name} {main_indicator}"
                    else:
                        main_tab_text = display_name

                    self.tabs.setTabText(index, main_tab_text)

                else:
                    # Tab has no subtabs - simple status update
                    status = self.tab_status[tab_key]
                    indicator = status_indicators.get(status, '')

                    if indicator:
                        tab_text = f"{display_name} {indicator}"
                    else:
                        tab_text = display_name

                    self.tabs.setTabText(index, tab_text)

    def set_tab_running(self, tab_name, subtab_name=None):
        """Set a tab or subtab status to running."""
        self.update_tab_status(tab_name, 'running', subtab_name)

    def set_tab_success(self, tab_name, subtab_name=None):
        """Set a tab or subtab status to success."""
        self.update_tab_status(tab_name, 'success', subtab_name)

    def set_tab_error(self, tab_name, subtab_name=None):
        """Set a tab or subtab status to error."""

        self.update_tab_status(tab_name, 'error', subtab_name)

    def set_tab_idle(self, tab_name, subtab_name=None):
        """Set a tab or subtab status to idle."""
        self.update_tab_status(tab_name, 'idle', subtab_name)

    def create_presimulation_tab(self):
        """
        Create the Presimulation tab with two subtabs for molecule generation and isomorphism computation.

        This method initializes the Presimulation tab in the GUI, which contains two subtabs:
        1. Generate Molfiles: Interface for generating molecular structure files
        2. Compute Isomorphism: Interface for computing structural similarities between molecules

        The tab also includes a console output area at the bottom for displaying execution logs
        and command outputs related to presimulation operations.

        Parameters:
            self: The POMSimulatorGUI instance

        Returns:
            None
        """
        tab = QWidget()
        self.tabs.addTab(tab, "Presimulation")
        layout = QVBoxLayout(tab)

        # Create subtabs
        subtabs = QTabWidget()
        layout.addWidget(subtabs)

        # Section 1: Generate Molfiles
        gen_mol_widget = QWidget()
        gen_mol_layout = QVBoxLayout(gen_mol_widget)
        self.create_gen_mol_section(gen_mol_layout)
        subtabs.addTab(gen_mol_widget, "Generate Molfiles")

        # Section 2: Compute Isomorphism
        iso_widget = QWidget()
        iso_layout = QVBoxLayout(iso_widget)
        self.create_iso_section(iso_layout)
        subtabs.addTab(iso_widget, "Compute Isomorphism")

    def create_gen_mol_section(self, layout):
        """
        Create the Generate Molfiles section of the Presimulation tab.

        This method builds the UI components for the Generate Molfiles section,
        which allows users to configure and run the molecule file generation process.
        It creates input fields for directory selection, output format options,
        and adds a run button to execute the generation script.

        Parameters:
            layout (QVBoxLayout): The parent layout where this section will be added.
                                 All UI components created by this method will be
                                 organized within this layout.

        Returns:
            None: This method modifies the provided layout in-place by adding
                  UI components to it.
        """

        # Parameters from generate_mol_file.py
        group = QGroupBox("Generate Molfiles Parameters")
        form_layout = QVBoxLayout()

        # Input directory
        input_layout = QHBoxLayout()
        input_layout.addWidget(QLabel("ADF Folder:"))
        self.adf_inputs_dir = QLineEdit()
        input_layout.addWidget(self.adf_inputs_dir)
        browse_btn = QPushButton("Browse")
        browse_btn.setToolTip("Browse to the directory where the input files are located.")
        browse_btn.clicked.connect(lambda: self.browse_directory(self.adf_inputs_dir))
        input_layout.addWidget(browse_btn)
        form_layout.addLayout(input_layout)

        self.register_ui_element(self.adf_inputs_dir, 'adf_inputs_dir')

        # Output directory
        output_layout = QHBoxLayout()
        output_layout.addWidget(QLabel("Molfile Directory:"))
        self.mol_folder = QLineEdit()
        output_layout.addWidget(self.mol_folder)
        browse_btn = QPushButton("Browse")
        browse_btn.setToolTip("Browse to the directory where the generated mol files will be saved.")
        browse_btn.clicked.connect(lambda: self.browse_directory(self.mol_folder))
        output_layout.addWidget(browse_btn)
        form_layout.addLayout(output_layout)

        self.register_ui_element(self.mol_folder, 'mol_folder')

        # # File format
        # format_layout = QHBoxLayout()
        # format_layout.addWidget(QLabel("Output Format:"))
        # self.gen_mol_format = QComboBox()
        # self.gen_mol_format.addItems(["MOL"])
        # format_layout.addWidget(self.gen_mol_format)
        # form_layout.addLayout(format_layout)

        group.setLayout(form_layout)
        layout.addWidget(group)

        # Run button
        run_btn = QPushButton("▶ Run Generate Molfiles")
        run_btn.setToolTip("Run the mol file generation script.")
        run_btn.setStyleSheet(ButtonStylesheets.RUN_BUTTON)

        self.gen_mol_run_btn = run_btn
        self.gen_mol_run_btn.clicked.connect(self.run_generate_mol_files)
        layout.addWidget(run_btn)

        # Stop button
        stop_button = QPushButton("⏹ Stop")
        stop_button.setToolTip("Stop the current operation.")
        stop_button.setStyleSheet(ButtonStylesheets.STOP_BUTTON)
        self.gen_mol_stop_btn = stop_button
        self.gen_mol_stop_btn.clicked.connect(self.stop_function)
        self.gen_mol_stop_btn.setEnabled(False)  # Initially disabled

        layout.addWidget(self.gen_mol_stop_btn)

    def run_generate_mol_files(self):
        """Run the generate mol files function"""
        self.clear_console()
        try:
            from utilities.generate_mol_file_gui import generate_molfile
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import generate_molfile: {str(e)}")
            self.set_tab_error('presimulation', 'generate_mol')
            return
        try:
            config_dict = self.get_gui_params()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to get GUI parameters: {str(e)}")
            self.set_tab_error('presimulation', 'generate_mol')
            return
        if not config_dict['Preparation']['adf_inputs_dir'] or not config_dict['Preparation']['mol_folder']:
            QMessageBox.critical(self, "Error", "Please select both ADF and Molfile directories.")
            self.set_tab_error('presimulation', 'generate_mol')
            return

        try:
            self.update_console("Started generate mol files function", "presimulation")
            self.set_tab_running('presimulation', 'generate_mol')

            # Create and start function runner
            self.mol_file_runner = POMSim_func_runner(generate_molfile, config_dict)
            self.mol_file_runner.stopped.connect(self.on_function_stopped)
            self.mol_file_runner.finished.connect(self.on_function_finished)
            self.mol_file_runner.error.connect(self.on_function_error)

            # Connect console output to update the console
            self.mol_file_runner.console_output.connect(
                lambda text: self.update_console(text, "presimulation"))
            self.gen_mol_run_btn.setEnabled(False)
            self.gen_mol_stop_btn.setEnabled(True)
            # Start the function
            self.mol_file_runner.start()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to start mol file generation: {str(e)}")
            self.set_tab_error('presimulation', 'generate_mol')
            self.gen_mol_run_btn.setEnabled(True)
            self.gen_mol_stop_btn.setEnabled(False)

    def create_iso_section(self, layout):
        """
        Create the Compute Isomorphism section of the Presimulation tab.

        This method builds the UI components for the Compute Isomorphism section,
        which allows users to configure and run the molecular isomorphism computation.
        It creates input fields for molecule directory selection, output file path,
        cores amount, and adds a run button to execute the
        isomorphism computation directly.

        Parameters:
            layout (QVBoxLayout): The parent layout where this section will be added.
                                 All UI components created by this method will be
                                 organized within this layout.

        Returns:
            None: This method modifies the provided layout in-place by adding
                  UI components to it.
        """
        # Parameters from compute_isomorphism.py
        group = QGroupBox("Compute Isomorphism Parameters")
        form_layout = QVBoxLayout()

        # System
        system_layout = QHBoxLayout()
        system_layout.addWidget(QLabel("System:"))
        self.POM_system = self.get_or_create_widget('POM_system', QLineEdit)
        self.POM_system.setText("")
        system_layout.addWidget(self.POM_system)
        form_layout.addLayout(system_layout)

        self.register_ui_element(self.POM_system, 'POM_system')

        # Molecule directory
        mol_layout = QHBoxLayout()
        mol_layout.addWidget(QLabel("Molfile Directory:"))
        self.mol_folder = QLineEdit()
        mol_layout.addWidget(self.mol_folder)
        browse_btn = QPushButton("Browse")
        browse_btn.setToolTip("Select the directory containing the molfiles")
        browse_btn.clicked.connect(lambda: self.browse_directory(self.mol_folder))
        mol_layout.addWidget(browse_btn)
        form_layout.addLayout(mol_layout)

        self.register_ui_element(self.mol_folder, 'mol_folder')

        # Output dir
        output_layout = QHBoxLayout()
        output_layout.addWidget(QLabel("Output Directory:"))
        self.output_path = QLineEdit()
        output_layout.addWidget(self.output_path)
        browse_btn = QPushButton("Browse")
        browse_btn.setToolTip("Select the directory to save the isomorphism matrix")
        browse_btn.clicked.connect(lambda: self.browse_directory(self.output_path))
        output_layout.addWidget(browse_btn)
        form_layout.addLayout(output_layout)

        self.register_ui_element(self.output_path, 'output_path')

        # Isomorphism cores
        cores_layout = QHBoxLayout()
        cores_layout.addWidget(QLabel("Isomorphism cores:"))
        self.iso_cores = CustomDoubleSpinBox()
        self.iso_cores.setRange(int(0), int(cpu_count() - 1))
        self.iso_cores.setValue(int(1))
        self.iso_cores.setSingleStep(int(1))
        cores_layout.addWidget(self.iso_cores)
        form_layout.addLayout(cores_layout)

        group.setLayout(form_layout)
        layout.addWidget(group)

        # Run button
        run_btn = QPushButton("▶ Run Compute Isomorphism")
        run_btn.setToolTip("Run the isomorphism calculation script.")
        run_btn.setStyleSheet(ButtonStylesheets.RUN_BUTTON)
        self.comp_iso_run_btn = run_btn
        self.comp_iso_run_btn.clicked.connect(self.run_compute_isomorphism)
        layout.addWidget(run_btn)

        # Stop button
        stop_button = QPushButton("⏹ Stop")
        stop_button.setToolTip("Stop the current operation.")
        stop_button.setStyleSheet(ButtonStylesheets.STOP_BUTTON)
        self.comp_iso_stop_btn = stop_button
        self.comp_iso_stop_btn.clicked.connect(self.stop_function)
        self.comp_iso_stop_btn.setEnabled(False)  # Initially disabled
        layout.addWidget(self.comp_iso_stop_btn)

    def run_compute_isomorphism(self):
        """Run the compute isomorphism function"""
        self.clear_console()
        try:
            from utilities.compute_isomorphism_gui import compute_isomorphism
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import compute isomorphism: {str(e)}")
            self.set_tab_error('presimulation', 'compute_iso')
            return
        try:
            config_dict = self.get_gui_params()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to get GUI parameters: {str(e)}")
            self.set_tab_error('presimulation', 'compute_iso')
            return
        if not config_dict['Preparation']['adf_inputs_dir'] or not config_dict['Preparation']['mol_folder'] or not \
                config_dict['Preparation']['output_path'] or not config_dict['Preparation']['POM_system']:
            QMessageBox.critical(self, "Error", "Please select both ADF and Molfile directories.")
            self.set_tab_error('presimulation', 'compute_iso')
            return

        try:
            self.update_console("Started compute isomorphism function", "presimulation")
            self.set_tab_running('presimulation', 'compute_iso')
            # Create and start function runner
            self.isomorph_runner = POMSim_func_runner(compute_isomorphism, config_dict)
            self.isomorph_runner.finished.connect(self.on_function_finished)
            self.isomorph_runner.error.connect(self.on_function_error)
            self.isomorph_runner.stopped.connect(self.on_function_stopped)
            # Connect console output to update the console
            self.isomorph_runner.console_output.connect(
                lambda text: self.update_console(text, "presimulation"))

            self.comp_iso_run_btn.setEnabled(False)
            self.comp_iso_stop_btn.setEnabled(True)
            # Start the function
            self.isomorph_runner.start()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to start compute isomorphism: {str(e)}")
            self.set_tab_error('presimulation', 'compute_iso')

    def create_simulation_tab(self):
        tab = QWidget()
        self.tabs.addTab(tab, "Simulation")
        layout = QVBoxLayout(tab)

        # Create subtabs
        subtabs = QTabWidget()
        layout.addWidget(subtabs)

        # Subtab compute lgkf
        comp_lgkf_widget = QWidget()
        comp_lgkf_layout = QVBoxLayout(comp_lgkf_widget)
        self.create_compute_lgkf_subtab(comp_lgkf_layout)
        subtabs.addTab(comp_lgkf_widget, "Compute lgkf")

        # Subtab generate CRN
        # gen_crn_widget = QWidget()
        # gen_crn_layout = QVBoxLayout(gen_crn_widget)
        # self.create_generate_crn_subtab(gen_crn_layout)
        # subtabs.addTab(gen_crn_widget, "Generate CRN")

    def create_compute_lgkf_subtab(self, layout):
        """
                Create the Simulation tab with IPA/HPA selection options.

                This method initializes the Simulation tab in the GUI, which allows users to
                choose between IPA (Isopolyanion) and HPA (Heteropolyanion) simulation types.
                The tab includes a dropdown for simulation type selection, a dynamic parameter
                section that updates based on the selected simulation type, and a console
                output area for displaying execution logs.

                The method creates the basic structure of the tab and calls update_sim_parameters()
                to populate the parameter section based on the initially selected simulation type.

                Parameters:
                    self: The POMSimulatorGUI instance

                Returns:
                    None: This method doesn't return a value but modifies the GUI by adding
                          a new tab to the main tab widget.
                """

        # Parameters container
        self.sim_params_container = QWidget()
        self.sim_params_layout = QVBoxLayout(self.sim_params_container)
        layout.addWidget(self.sim_params_container)

        # Initialize parameters
        self.update_sim_parameters()

        # Add console at bottom

    def _build_sim_io_group(self, layout):
        """Build the shared Input/Output File Managing group for the simulation tab.

        Creates System, ADF Folder, Mol Folder, and Output Path fields and registers
        them as shared UI elements. Adds the resulting QGroupBox to *layout*.
        """
        system_group = QGroupBox("Input Output File Managing")
        system_layout = QVBoxLayout()

        # System
        system_row = QHBoxLayout()
        system_row.addWidget(QLabel("System:"))
        self.POM_system = QLineEdit()
        system_row.addWidget(self.POM_system)
        system_layout.addLayout(system_row)
        self.register_ui_element(self.POM_system, "POM_system")

        # ADF Folder
        adf_row = QHBoxLayout()
        adf_row.addWidget(QLabel("ADF Folder:"))
        self.adf_inputs_dir = QLineEdit()
        adf_row.addWidget(self.adf_inputs_dir)
        browse_adf = QPushButton("Browse")
        browse_adf.clicked.connect(lambda: self.browse_directory(self.adf_inputs_dir))
        adf_row.addWidget(browse_adf)
        system_layout.addLayout(adf_row)
        self.register_ui_element(self.adf_inputs_dir, "adf_inputs_dir")

        # Mol Folder
        mol_row = QHBoxLayout()
        mol_row.addWidget(QLabel("Mol Folder:"))
        self.mol_folder = QLineEdit()
        mol_row.addWidget(self.mol_folder)
        browse_mol = QPushButton("Browse")
        browse_mol.clicked.connect(lambda: self.browse_directory(self.mol_folder))
        mol_row.addWidget(browse_mol)
        system_layout.addLayout(mol_row)
        self.register_ui_element(self.mol_folder, "mol_folder")

        # Output Path
        output_row = QHBoxLayout()
        output_row.addWidget(QLabel("Output Path:"))
        self.output_path = QLineEdit()
        output_row.addWidget(self.output_path)
        browse_output = QPushButton("Browse")
        browse_output.clicked.connect(lambda: self.browse_directory(self.output_path))
        output_row.addWidget(browse_output)
        system_layout.addLayout(output_row)
        self.register_ui_element(self.output_path, "output_path")

        system_group.setLayout(system_layout)
        layout.addWidget(system_group)

    def _build_sim_operation_group(self, layout):
        """Build the shared Operation Parameters group for the simulation tab.

        Creates Cores, Batch Size, Sample Percentage, and Sample Type fields.
        Adds the resulting QGroupBox to *layout*.
        """
        sim_settings_group = QGroupBox("Operation Parameters")
        sim_settings_layout = QVBoxLayout()

        # Cores
        cores_row = QHBoxLayout()
        cores_row.addWidget(QLabel("Cores:"))
        self.sim_cores = QSpinBox()
        self.sim_cores.setRange(1, cpu_count() - 1)
        self.sim_cores.setValue(1)
        cores_row.addWidget(self.sim_cores)
        sim_settings_layout.addLayout(cores_row)

        # Batch Size
        batch_row = QHBoxLayout()
        batch_row.addWidget(QLabel("Batch Size:"))
        self.sim_batch_size = QSpinBox()
        self.sim_batch_size.setRange(1, 1000)
        self.sim_batch_size.setValue(200)
        batch_row.addWidget(self.sim_batch_size)
        sim_settings_layout.addLayout(batch_row)

        # Sample Percentage
        sample_perc_row = QHBoxLayout()
        sample_perc_row.addWidget(QLabel("Sample Percentage:"))
        self.sim_sample_perc = CustomDoubleSpinBox()
        self.sim_sample_perc.setRange(0, 100)
        self.sim_sample_perc.setDecimals(8)
        self.sim_sample_perc.setSingleStep(0.1)
        self.sim_sample_perc.setValue(10)
        sample_perc_row.addWidget(self.sim_sample_perc)
        sim_settings_layout.addLayout(sample_perc_row)

        # Sample Type
        sample_type_row = QHBoxLayout()
        sample_type_row.addWidget(QLabel("Sample Type:"))
        self.sim_sample_type = QComboBox()
        self.sim_sample_type.addItems(["random", "all"])
        sample_type_row.addWidget(self.sim_sample_type)
        sim_settings_layout.addLayout(sample_type_row)

        sim_settings_group.setLayout(sim_settings_layout)
        layout.addWidget(sim_settings_group)

    def _build_sim_common_chemical_group(self, layout):
        """Build the shared Chemical Parameters section common to both IPA and HPA.

        Creates Use Isomorphism, Energy Threshold, Proton Number, Reference Types
        grid, Ionic Strength, Temperature, pH Range, and MSCE Solver fields.
        Returns the QVBoxLayout of the Chemical parameters group so callers can
        append mode-specific concentration fields before finalising the group.
        """
        adv_settings_group = QGroupBox("Chemical parameters")
        adv_settings_layout = QVBoxLayout()

        # Use Isomorphism
        iso_row = QHBoxLayout()
        self.use_isomorphism = QCheckBox("Use Isomorphism")
        self.use_isomorphism.setChecked(True)
        iso_row.addWidget(self.use_isomorphism)
        adv_settings_layout.addLayout(iso_row)

        # Energy Threshold
        energy_row = QHBoxLayout()
        energy_row.addWidget(QLabel("Energy Threshold:"))
        self.energy_threshold = CustomDoubleSpinBox()
        self.energy_threshold.setRange(-500, 500)
        self.energy_threshold.setDecimals(2)
        self.energy_threshold.setSingleStep(0.1)
        self.energy_threshold.setValue(20)
        energy_row.addWidget(self.energy_threshold)
        adv_settings_layout.addLayout(energy_row)

        # Proton Number
        proton_row = QHBoxLayout()
        proton_row.addWidget(QLabel("Proton Number:"))
        self.proton_numb = QSpinBox()
        self.proton_numb.setRange(0, 100)
        self.proton_numb.setValue(0)
        proton_row.addWidget(self.proton_numb)
        adv_settings_layout.addLayout(proton_row)

        # Reference Types
        ref_types_group = QGroupBox("Reference Types")
        ref_types_layout = QGridLayout()
        self.ref_types = {}
        num_columns = int(len(reaction_references) // 4) + 1
        for ii, ref_type in enumerate(reaction_references):
            row = ii // num_columns
            col = ii % num_columns
            checkbox = QCheckBox(ref_type)
            checkbox.setChecked(True)
            self.ref_types[ref_type] = checkbox
            ref_types_layout.addWidget(checkbox, row, col)
        ref_types_group.setLayout(ref_types_layout)
        adv_settings_layout.addWidget(ref_types_group)

        # Ionic Strength
        i_row = QHBoxLayout()
        i_row.addWidget(QLabel("Ionic Strength (I):"))
        self.i_s = CustomDoubleSpinBox()
        self.i_s.setRange(0.001, 1)
        self.i_s.setDecimals(4)
        self.i_s.setSingleStep(0.001)
        self.i_s.setValue(0.1)
        i_row.addWidget(self.i_s)
        adv_settings_layout.addLayout(i_row)

        # Temperature
        temp_row = QHBoxLayout()
        temp_row.addWidget(QLabel("Temperature (K):"))
        self.temp = CustomDoubleSpinBox()
        self.temp.setRange(0, 600)
        self.temp.setValue(298.15)
        temp_row.addWidget(self.temp)
        adv_settings_layout.addLayout(temp_row)

        # pH Range
        ph_group = QGroupBox("pH Range")
        ph_layout = QVBoxLayout()

        min_ph_row = QHBoxLayout()
        min_ph_row.addWidget(QLabel("Min pH:"))
        self.sim_min_ph = CustomDoubleSpinBox()
        self.sim_min_ph.setRange(-10, 100)
        self.sim_min_ph.setDecimals(2)
        self.sim_min_ph.setSingleStep(0.01)
        self.sim_min_ph.setValue(0)
        min_ph_row.addWidget(self.sim_min_ph)
        ph_layout.addLayout(min_ph_row)

        max_ph_row = QHBoxLayout()
        max_ph_row.addWidget(QLabel("Max pH:"))
        self.sim_max_ph = CustomDoubleSpinBox()
        self.sim_max_ph.setRange(-10, 100)
        self.sim_max_ph.setDecimals(2)
        self.sim_max_ph.setSingleStep(0.01)
        self.sim_max_ph.setValue(35)
        max_ph_row.addWidget(self.sim_max_ph)
        ph_layout.addLayout(max_ph_row)

        step_ph_row = QHBoxLayout()
        step_ph_row.addWidget(QLabel("Step pH:"))
        self.sim_step_ph = CustomDoubleSpinBox()
        self.sim_step_ph.setRange(0.01, 1)
        self.sim_step_ph.setDecimals(2)
        self.sim_step_ph.setSingleStep(0.01)
        self.sim_step_ph.setValue(0.5)
        step_ph_row.addWidget(self.sim_step_ph)
        ph_layout.addLayout(step_ph_row)

        ph_group.setLayout(ph_layout)
        adv_settings_layout.addWidget(ph_group)


        # Return the group and its layout so callers can append mode-specific fields
        return adv_settings_group, adv_settings_layout

    def _build_sim_internal_group(self, layout):
        """Build the shared Internal Parameters group for the simulation tab.

        Creates Restrain Addition, Restrain Condensation, Include Dimerization,
        Force Stoich, and Adjust Hydration Protons fields.
        Adds the resulting QGroupBox to *layout*.
        """
        int_settings_group = QGroupBox("Internal parameters")
        int_settings_layout = QVBoxLayout()

        # Restrain Addition
        rest_add_row = QHBoxLayout()
        rest_add_row.addWidget(QLabel("Restrain Addition:"))
        self.rest_add = CustomDoubleSpinBox()
        self.rest_add.setRange(0, 200)
        self.rest_add.setValue(0)
        rest_add_row.addWidget(self.rest_add)
        int_settings_layout.addLayout(rest_add_row)

        # Restrain Condensation
        rest_cond_row = QHBoxLayout()
        rest_cond_row.addWidget(QLabel("Restrain Condensation:"))
        self.rest_cond = CustomDoubleSpinBox()
        self.rest_cond.setRange(0, 200)
        self.rest_cond.setValue(0)
        rest_cond_row.addWidget(self.rest_cond)
        int_settings_layout.addLayout(rest_cond_row)

        # Include dimerization
        inc_dim_row = QHBoxLayout()
        self.include_dimerization = QCheckBox("Include dimerization")
        self.include_dimerization.setChecked(True)
        inc_dim_row.addWidget(self.include_dimerization)
        int_settings_layout.addLayout(inc_dim_row)

        # Force stoich
        force_sto_row = QHBoxLayout()
        force_sto_row.addWidget(QLabel("Force stoich (comma-separated integers):"))
        self.force_sto = QLineEdit()
        self.force_sto.setPlaceholderText("e.g., 1,2,3")
        self.force_sto.setToolTip("Enter a comma-separated list of integers")
        force_sto_row.addWidget(self.force_sto)
        int_settings_layout.addLayout(force_sto_row)

        # Adjust protons hydration
        adj_prot_row = QHBoxLayout()
        self.adj_prot = QCheckBox("Adjust hydration protons")
        self.adj_prot.setChecked(True)
        adj_prot_row.addWidget(self.adj_prot)
        int_settings_layout.addLayout(adj_prot_row)

        int_settings_group.setLayout(int_settings_layout)
        layout.addWidget(int_settings_group)

    def update_sim_parameters(self):
        """
        Update the simulation parameters UI based on the selected simulation type.
        """
        # Clear existing parameters
        while self.sim_params_layout.count():
            child = self.sim_params_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
            elif child.layout():
                # Clear child layouts too
                while child.layout().count():
                    subchild = child.layout().takeAt(0)
                    if subchild.widget():
                        subchild.widget().deleteLater()

        # Create new parameters based on type
        pom_type = getattr(self, 'global_mode', 'IPA')

        # Create a scroll area to contain all parameters
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)

        # --- Shared groups (identical for IPA and HPA) ---
        self._build_sim_io_group(scroll_layout)
        self._build_sim_operation_group(scroll_layout)

        # --- Chemical parameters group (common base + mode-specific concentration) ---
        adv_settings_group, adv_settings_layout = self._build_sim_common_chemical_group(scroll_layout)

        if pom_type == "IPA":
            # IPA-specific: single initial concentration + single reference compound
            c0_row = QHBoxLayout()
            c0_row.addWidget(QLabel("Initial Concentration (C0):"))
            self.sim_c0 = CustomDoubleSpinBox()
            self.sim_c0.setRange(0.0000001, 5.0)
            self.sim_c0.setDecimals(8)
            self.sim_c0.setSingleStep(0.001)
            self.sim_c0.setValue(0.01)
            c0_row.addWidget(self.sim_c0)
            adv_settings_layout.addLayout(c0_row)

            ref_comp_row = QHBoxLayout()
            ref_comp_row.addWidget(QLabel("Reference Compound:"))
            self.ref_compound = QLineEdit()
            self.ref_compound.setText("")
            ref_comp_row.addWidget(self.ref_compound)
            adv_settings_layout.addLayout(ref_comp_row)

        elif pom_type == "HPA":
            # HPA-specific: metal + heteroatom concentrations and reference compounds
            CM_row = QHBoxLayout()
            CM_row.addWidget(QLabel("Initial Metal Concentration (CM):"))
            self.sim_CM = CustomDoubleSpinBox()
            self.sim_CM.setRange(0.0000001, 5.0)
            self.sim_CM.setDecimals(8)
            self.sim_CM.setSingleStep(0.0001)
            self.sim_CM.setValue(0.005)
            CM_row.addWidget(self.sim_CM)
            adv_settings_layout.addLayout(CM_row)

            CX_row = QHBoxLayout()
            CX_row.addWidget(QLabel("Initial Heteroatom Concentration (CX):"))
            self.sim_CX = CustomDoubleSpinBox()
            self.sim_CX.setRange(0.0000001, 5.0)
            self.sim_CX.setDecimals(8)
            self.sim_CX.setSingleStep(0.0001)
            self.sim_CX.setValue(0.005)
            CX_row.addWidget(self.sim_CX)
            adv_settings_layout.addLayout(CX_row)

            refM_comp_row = QHBoxLayout()
            refM_comp_row.addWidget(QLabel("Reference Metal Compound:"))
            self.sim_M_ref_compound = QLineEdit()
            self.sim_M_ref_compound.setText("")
            refM_comp_row.addWidget(self.sim_M_ref_compound)
            adv_settings_layout.addLayout(refM_comp_row)

            refX_comp_row = QHBoxLayout()
            refX_comp_row.addWidget(QLabel("Reference Heteroatom Compound:"))
            self.sim_X_ref_compound = QLineEdit()
            self.sim_X_ref_compound.setText("")
            refX_comp_row.addWidget(self.sim_X_ref_compound)
            adv_settings_layout.addLayout(refX_comp_row)

        adv_settings_group.setLayout(adv_settings_layout)
        scroll_layout.addWidget(adv_settings_group)

        # --- Internal parameters (identical for IPA and HPA) ---
        self._build_sim_internal_group(scroll_layout)

        # Finalise scroll area
        scroll_widget.setLayout(scroll_layout)
        scroll_area.setWidget(scroll_widget)
        self.sim_params_layout.addWidget(scroll_area)

        # Run / Stop buttons (mode-specific attribute names for backward compatibility)
        run_btn = QPushButton(f"▶ Run {pom_type} Simulation")
        run_btn.setToolTip("Run the simulation.")
        run_btn.setStyleSheet(ButtonStylesheets.RUN_BUTTON)
        run_btn.clicked.connect(self.run_compute_lgkf)

        stop_btn = QPushButton("⏹ Stop")
        stop_btn.setToolTip("Stop the current simulation.")
        stop_btn.setStyleSheet(ButtonStylesheets.STOP_BUTTON)
        stop_btn.setEnabled(False)
        stop_btn.clicked.connect(self.stop_function)

        if pom_type == "IPA":
            self.ipa_comp_lgkf_run_btn = run_btn
            self.ipa_comp_lgkf_stop_btn = stop_btn
        else:
            self.hpa_comp_lgkf_run_btn = run_btn
            self.hpa_comp_lgkf_stop_btn = stop_btn

        self.sim_params_layout.addWidget(run_btn)
        self.sim_params_layout.addWidget(stop_btn)

    def run_compute_lgkf(self):

        """Run the compute_lgkf function"""
        self.clear_console()

        try:
            if self.pom_type == "IPA":
                from simulations.simulation_gui import simulation_ipa
                simulation_function = simulation_ipa
                tab_context = "IPA Simulation"
            elif self.pom_type == "HPA":
                from simulations.simulation_gui import simulation_hpa
                simulation_function = simulation_hpa
                tab_context = "HPA Simulation"
            else:
                QMessageBox.critical(self, "Error", "Invalid simulation type selected.")
                self.set_tab_error('simulation', 'compute_lgkf')
                return
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import simulation function: {str(e)}")
            self.set_tab_error('simulation', 'compute_lgkf')
            return

        try:
            config_dict = self.get_gui_params()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to get GUI parameters: {str(e)}")
            self.set_tab_error('simulation', 'compute_lgkf')
            return

        if not config_dict['Preparation']['output_path'] or not config_dict['Preparation']['mol_folder'] or not \
                config_dict['Preparation']['POM_system'] or not config_dict['Preparation']['adf_inputs_dir'] \
                or not config_dict['Simulation']['ref_compound']:
            QMessageBox.critical(self, "Error", "No ADF and Molfile directories or POM system or output path.")
            self.set_tab_error('simulation', 'compute_lgkf')
            return

        try:
            self.update_console("Started compute_lgkf function", "simulation")
            self.set_tab_running('simulation', 'compute_lgkf')

            # Create and start function runner
            self.simulation_runner = POMSim_func_runner(simulation_function, config_dict)
            self.simulation_runner.finished.connect(self.on_function_finished)
            self.simulation_runner.error.connect(self.on_function_error)
            self.simulation_runner.stopped.connect(self.on_function_stopped)
            # Connect console output to update the console
            self.simulation_runner.console_output.connect(
                lambda text: self.update_console(text, "simulation"))
            if self.pom_type == "IPA":
                self.ipa_comp_lgkf_run_btn.setEnabled(False)
                self.ipa_comp_lgkf_stop_btn.setEnabled(True)
            elif self.pom_type == "HPA":
                self.hpa_comp_lgkf_run_btn.setEnabled(False)
                self.hpa_comp_lgkf_stop_btn.setEnabled(True)
            # Start the function
            self.simulation_runner.start()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to start {tab_context}: {str(e)}")
            self.set_tab_error('simulation', 'compute_lgkf')
            # Reset button states on error
            if self.pom_type == "IPA":
                self.ipa_comp_lgkf_run_btn.setEnabled(True)
                self.ipa_comp_lgkf_stop_btn.setEnabled(False)
            elif self.pom_type == "HPA":
                self.hpa_comp_lgkf_run_btn.setEnabled(True)
                self.hpa_comp_lgkf_stop_btn.setEnabled(False)

    def create_generate_crn_subtab(self, layout):
        """
        Create the Generate CRN tab with input and output file selection options.

        This method initializes the Generate CRN tab in the GUI, which allows users to
        select input and output files for CRN generation. The tab includes a file"""

        # todo: Implement CRN generation logic here

    def create_scaling_tab(self):
        """
        Create the Scaling tab in the GUI.

        This method initializes the Scaling tab, which provides an interface for users
        to configure and run the scaling constants functionality. The tab includes
        input fields for selecting files to process, parameter configuration options,
        a run button to execute the scaling script, and a console output area for
        displaying execution logs.

        The scaling functionality allows users to adjust constants used in the simulation
        based on input data, which is critical for calibrating the simulation to match
        experimental results.

        Parameters:
            self: The POMSimulatorGUI instance containing the tab widget where
                  this tab will be added.

        Returns:
            None: This method doesn't return a value but modifies the GUI by adding
                  a new tab to the main tab widget.
        """
        tab = QWidget()
        self.tabs.addTab(tab, "Scaling")
        layout = QVBoxLayout(tab)

        # Add scaling parameters
        group = QGroupBox("Scaling Parameters")
        form_layout = QVBoxLayout()

        # Parameters from scale_constants.py
        pom_system_layout = QHBoxLayout()
        pom_system_layout.addWidget(QLabel("POM System:"))
        self.POM_system = QLineEdit()
        pom_system_layout.addWidget(self.POM_system)
        form_layout.addLayout(pom_system_layout)

        self.register_ui_element(self.POM_system, "POM_system")

        # Output path layout
        output_path_layout = QHBoxLayout()
        output_path_layout.addWidget(QLabel("Output Path:"))
        self.output_path = QLineEdit()
        output_path_layout.addWidget(self.output_path)
        output_browse_btn = QPushButton("Browse")
        output_browse_btn.setToolTip("Browse to the directory where the input files are located.")
        output_browse_btn.clicked.connect(lambda: self.browse_directory(self.output_path))
        output_path_layout.addWidget(output_browse_btn)
        form_layout.addLayout(output_path_layout)

        self.register_ui_element(self.output_path, "output_path")

        # Scaling mode

        scaling_mode_layout = QHBoxLayout()
        scaling_mode_layout.addWidget(QLabel("Scaling Mode:"))
        self.scaling_mode = QComboBox()
        self.scaling_mode.addItems(allowed_scaling_modes)
        scaling_mode_layout.addWidget(self.scaling_mode)
        form_layout.addLayout(scaling_mode_layout)

        # Experimental set layout
        exp_set_layout = QHBoxLayout()
        exp_set_layout.addWidget(QLabel("Experimental Set:"))
        self.exp_set = QComboBox()
        self.exp_set.addItems(list(experimental_constants.keys()))
        exp_set_layout.addWidget(self.exp_set)
        form_layout.addLayout(exp_set_layout)

        # Add a display area for the constants
        constants_group = QGroupBox("Experimental Constants")
        constants_layout = QVBoxLayout()

        # Create a text area to display the constants
        self.constants_display = QTextEdit()
        self.constants_display.setReadOnly(True)
        self.constants_display.setMinimumHeight(40)

        # Apply initial theme-aware styling
        if self.dark_theme:
            self.constants_display.setStyleSheet("""
                QTextEdit {
                    background-color: #2b2b2b;
                    color: #ffffff;
                    font-family: monospace;
                    border: 1px solid #555555;
                }
            """)
        else:
            self.constants_display.setStyleSheet("""
                QTextEdit {
                    background-color: #f8f8f8;
                    color: #000000;
                    font-family: monospace;
                    border: 1px solid #cccccc;
                }
            """)

        constants_layout.addWidget(self.constants_display)

        # Update the display when a different set is selected
        self.exp_set.currentIndexChanged.connect(self.update_constants_display)
        constants_group.setLayout(constants_layout)
        form_layout.addWidget(constants_group)

        # Initial update of the display with the first set
        self.update_constants_display()

        # Additional scaling parameters...

        group.setLayout(form_layout)
        layout.addWidget(group)

        # Run button
        run_btn = QPushButton("▶ Run Scaling")
        run_btn.setToolTip("Run the scaling.")
        run_btn.setStyleSheet(ButtonStylesheets.RUN_BUTTON)
        run_btn.clicked.connect(self.run_scale_constants)
        self.scale_run_btn = run_btn
        layout.addWidget(self.scale_run_btn)

        # Stop button
        stop_button = QPushButton("⏹ Stop")
        stop_button.setToolTip("Stop the current operation.")
        stop_button.setStyleSheet(ButtonStylesheets.STOP_BUTTON)
        stop_button.setEnabled(False)  # Initially disabled
        self.scale_stop_btn = stop_button
        self.scale_stop_btn.clicked.connect(self.stop_function)
        layout.addWidget(self.scale_stop_btn)

        # Add console at bottom

    def update_constants_display(self):
        """
        Update the display of experimental constants when a different set is selected.
        """
        selected_set = self.exp_set.currentText()
        if selected_set in experimental_constants:
            constants = experimental_constants[selected_set]

            # Format the constants for display
            display_text = f"Constants for {selected_set}:\n\n"

            # Format each constant with its name and value
            for name, value in constants.items():
                if isinstance(value, (int, float)):
                    # Format numbers with appropriate precision
                    if isinstance(value, int):
                        formatted_value = str(value)
                    else:
                        formatted_value = f"{value:.4f}"
                    display_text += f"{name}: {formatted_value}\n"
                else:
                    display_text += f"{name}: {value}\n"

            self.constants_display.setText(display_text)
        else:
            self.constants_display.setText(f"No constants found for {selected_set}")

        # Apply theme-aware styling (same for both success and error cases)
        if self.dark_theme:
            self.constants_display.setStyleSheet("""
                QTextEdit {
                    background-color: #2b2b2b;
                    color: #ffffff;
                    font-family: monospace;
                    border: 1px solid #555555;
                }
            """)
        else:
            self.constants_display.setStyleSheet("""
                QTextEdit {
                    background-color: #f8f8f8;
                    color: #000000;
                    font-family: monospace;
                    border: 1px solid #cccccc;
                }
            """)

    def run_scale_constants(self):
        """Run the scaling of pomsimulator formation constants."""
        self.clear_console()
        try:
            from utilities.scale_constants import scale_constants
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import scale_constants. Error: {str(e)}")
            self.set_tab_error('scaling')

        try:
            config_dict = self.get_gui_params()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to retrieve GUI parameters. Error: {str(e)}")
            self.set_tab_error('scaling')
            return

        if not config_dict['Preparation']['output_path'] or not config_dict['Preparation']['POM_system']:
            QMessageBox.critical(self, "Error", "Please select an output path and a POM system.")
            self.set_tab_error('scaling')
            return

        try:
            self.update_console("Started scaling of pomsimulator formation constants...", "Scaling")
            self.set_tab_running('scaling')

            # Create and start function runner
            self.scaling_runner = POMSim_func_runner(scale_constants, config_dict)
            self.scaling_runner.finished.connect(self.on_function_finished)
            self.scaling_runner.error.connect(self.on_function_error)
            self.scaling_runner.stopped.connect(self.on_function_stopped)

            self.scaling_runner.console_output.connect(
                lambda text: self.update_console(text, "Scaling"))

            self.scale_run_btn.setEnabled(False)  # Disable run button initially
            self.scale_stop_btn.setEnabled(True)  # Enable stop button initially

            self.scaling_runner.start()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to run scaling. Error: {str(e)}")
            self.set_tab_error('Scaling')

            self.scale_run_btn.setEnabled(True)  # Enable run button
            self.scale_stop_btn.setEnabled(False)  # Disable stop button

    def create_speciation_tab(self):
        """
        Create the Speciation tab with two subtabs for different speciation analyses.

        This method initializes the Speciation tab in the GUI, which contains two subtabs:
        1. SpeciationA: For analyzing species distribution as a function of pH at a fixed concentration
        2. SpeciationN: For analyzing species distribution with varying parameters based on model index

        The tab includes a tabbed interface for the two speciation types and a console output
        area at the bottom for displaying execution logs related to speciation operations.

        Parameters:
            self: The POMSimulatorGUI instance containing the tab widget where
                  this tab will be added.

        Returns:
            None: This method doesn't return a value but modifies the GUI by adding
                  a new tab to the main tab widget.
        """
        tab = QWidget()
        self.tabs.addTab(tab, "Speciation")
        layout = QVBoxLayout(tab)

        # Create subtabs
        subtabs = QTabWidget()
        layout.addWidget(subtabs)

        # Create the speciation diagram tab
        speciation_diagram_widget = QWidget()
        speciation_diagram_layout = QVBoxLayout(speciation_diagram_widget)
        self.create_speciation_diagram_subtab(speciation_diagram_layout)
        subtabs.addTab(speciation_diagram_widget, "Speciation Diagram")

        # Create the speciation diagram tab
        phase_diagram_widget = QWidget()
        phase_diagram_layout = QVBoxLayout(phase_diagram_widget)
        self.create_phase_diagram_subtab(phase_diagram_layout)
        subtabs.addTab(phase_diagram_widget, "Phase Diagram")

    def create_speciation_diagram_subtab(self, layout):
        """
        Create the Speciation Diagram subtab with parameters for generating speciation diagrams.

        Args:
            layout (QVBoxLayout): The parent layout to add components to
        """
        # Parameters container
        self.speciation_params_container = QWidget()
        self.speciation_params_layout = QVBoxLayout(self.speciation_params_container)
        layout.addWidget(self.speciation_params_container)

        # Initialize parameters for Speciation
        self.update_speciation_parameters()

        # Add console at bottom

    def _build_spec_io_group(self, layout):
        """Build the shared Input/Output File Management group for speciation/phase tabs.

        Creates System and Output Path fields and registers them as shared UI elements.
        Adds the resulting QGroupBox to *layout*.
        """
        system_group = QGroupBox("Input/Output File Management")
        system_layout = QVBoxLayout()

        system_row = QHBoxLayout()
        system_row.addWidget(QLabel("System:"))
        self.POM_system = QLineEdit()
        system_row.addWidget(self.POM_system)
        system_layout.addLayout(system_row)
        self.register_ui_element(self.POM_system, "POM_system")

        output_row = QHBoxLayout()
        output_row.addWidget(QLabel("Output Path:"))
        self.output_path = QLineEdit()
        output_row.addWidget(self.output_path)
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(lambda: self.browse_directory(self.output_path))
        output_row.addWidget(browse_btn)
        system_layout.addLayout(output_row)
        self.register_ui_element(self.output_path, "output_path")

        system_group.setLayout(system_layout)
        layout.addWidget(system_group)

    def _build_spec_ph_group(self, ph_attr_prefix="spec"):
        """Build a pH Range QGroupBox with min/max/step spinboxes.

        Stores the created widgets as ``self.<ph_attr_prefix>_min_ph``,
        ``self.<ph_attr_prefix>_max_ph``, and ``self.<ph_attr_prefix>_step_ph``.

        Returns the configured QGroupBox.
        """
        ph_group = QGroupBox("pH Range")
        ph_layout = QVBoxLayout()

        min_ph_row = QHBoxLayout()
        min_ph_row.addWidget(QLabel("Min pH:"))
        min_ph = CustomDoubleSpinBox()
        min_ph.setRange(-10, 100)
        min_ph.setValue(0)
        min_ph_row.addWidget(min_ph)
        ph_layout.addLayout(min_ph_row)
        setattr(self, f"{ph_attr_prefix}_min_ph", min_ph)

        max_ph_row = QHBoxLayout()
        max_ph_row.addWidget(QLabel("Max pH:"))
        max_ph = CustomDoubleSpinBox()
        max_ph.setRange(-10, 100)
        max_ph.setValue(14)
        max_ph_row.addWidget(max_ph)
        ph_layout.addLayout(max_ph_row)
        setattr(self, f"{ph_attr_prefix}_max_ph", max_ph)

        step_ph_row = QHBoxLayout()
        step_ph_row.addWidget(QLabel("step pH:"))
        step_ph = CustomDoubleSpinBox()
        step_ph.setRange(0.01, 1)
        step_ph.setValue(0.5)
        step_ph.setDecimals(2)
        step_ph.setSingleStep(0.1)
        step_ph_row.addWidget(step_ph)
        ph_layout.addLayout(step_ph_row)
        setattr(self, f"{ph_attr_prefix}_step_ph", step_ph)

        ph_group.setLayout(ph_layout)
        return ph_group

    def _build_spec_labels_group(self, layout, labels_attr="spec_labels"):
        """Build the Species to Plot group with Labels File and Selected Labels fields.

        Stores the labels file widget as ``self.labels_file`` and the selected labels
        widget as ``self.<labels_attr>``. Registers both as shared UI elements.
        Adds the resulting QGroupBox to *layout*.
        """
        plot_list_group = QGroupBox("Species to Plot")
        plot_list_layout = QVBoxLayout()
        plot_list_group.setLayout(plot_list_layout)

        labels_file_row = QHBoxLayout()
        labels_file_row.addWidget(QLabel("Labels File:"))
        self.labels_file = QLineEdit()
        labels_file_row.addWidget(self.labels_file)
        labels_file_browse_btn = QPushButton("Browse")
        labels_file_browse_btn.clicked.connect(lambda: self.browse_file(self.labels_file))
        self.register_ui_element(self.labels_file, "labels_file")
        labels_file_row.addWidget(labels_file_browse_btn)
        plot_list_layout.addLayout(labels_file_row)

        selected_labels_row = QHBoxLayout()
        selected_labels_row.addWidget(QLabel("Selected Labels:"))
        selected_labels = QLineEdit()
        selected_labels.setText("all")
        selected_labels.setReadOnly(False)
        selected_labels_row.addWidget(selected_labels)
        setattr(self, labels_attr, selected_labels)
        select_labels_btn = QPushButton("Select Labels")
        select_labels_btn.clicked.connect(
            lambda: self.load_and_select_labels(self.labels_file, getattr(self, labels_attr))
        )
        self.register_ui_element(selected_labels, labels_attr)
        selected_labels_row.addWidget(select_labels_btn)
        plot_list_layout.addLayout(selected_labels_row)

        layout.addWidget(plot_list_group)

    def _build_spec_operation_group(self, layout):
        """Build the Operation Parameters group with Cores and Batch Size spinboxes.

        Stores widgets as ``self.spec_cores`` and ``self.spec_batch_size`` and
        registers them as shared UI elements. Adds the QGroupBox to *layout*.
        """
        params_group = QGroupBox("Operation Parameters")
        params_layout = QGridLayout()

        params_layout.addWidget(QLabel("Cores:"), 1, 0)
        self.spec_cores = QSpinBox()
        self.spec_cores.setRange(1, cpu_count() - 1)
        self.spec_cores.setValue(4)
        params_layout.addWidget(self.spec_cores, 1, 1)
        self.register_ui_element(self.spec_cores, "spec_cores")

        params_layout.addWidget(QLabel("Batch Size:"), 2, 0)
        self.spec_batch_size = QSpinBox()
        self.spec_batch_size.setRange(1, 1000)
        self.spec_batch_size.setValue(100)
        params_layout.addWidget(self.spec_batch_size, 2, 1)
        self.register_ui_element(self.spec_batch_size, "spec_batch_size")

        params_group.setLayout(params_layout)
        layout.addWidget(params_group)

    def update_speciation_parameters(self):
        """
        Update the speciation parameters UI based on the selected simulation type.
        """
        # Clear all existing widgets from the layout
        while self.speciation_params_layout.count():
            child = self.speciation_params_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

        pom_type = getattr(self, "global_mode", "IPA")

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)

        # --- Shared: I/O group ---
        self._build_spec_io_group(scroll_layout)

        # --- Chemical Parameters (mode-specific concentration + shared pH range) ---
        chem_group = QGroupBox("Chemical Parameters")
        chem_layout = QVBoxLayout()

        if pom_type == "IPA":
            C_row = QHBoxLayout()
            C_row.addWidget(QLabel("Initial Concentration (mol/L):"))
            self.spec_C0 = CustomDoubleSpinBox()
            self.spec_C0.setDecimals(6)
            self.spec_C0.setRange(0.000001, 10)
            self.spec_C0.setValue(0.1)
            self.spec_C0.setSingleStep(0.1)
            C_row.addWidget(self.spec_C0)
            chem_layout.addLayout(C_row)
        else:  # HPA
            CM_row = QHBoxLayout()
            CM_row.addWidget(QLabel("Initial Metal Concentration (mol/L):"))
            self.spec_CM = CustomDoubleSpinBox()
            self.spec_CM.setDecimals(6)
            self.spec_CM.setRange(0.000001, 10)
            self.spec_CM.setValue(0.1)
            self.spec_CM.setSingleStep(0.1)
            CM_row.addWidget(self.spec_CM)
            chem_layout.addLayout(CM_row)

            CX_row = QHBoxLayout()
            CX_row.addWidget(QLabel("Initial Heteroatom Concentration (mol/L):"))
            self.spec_CX = CustomDoubleSpinBox()
            self.spec_CX.setDecimals(6)
            self.spec_CX.setRange(0.000001, 10)
            self.spec_CX.setValue(0.1)
            self.spec_CX.setSingleStep(0.1)
            CX_row.addWidget(self.spec_CX)
            chem_layout.addLayout(CX_row)

        chem_layout.addWidget(self._build_spec_ph_group("spec"))
        chem_group.setLayout(chem_layout)
        scroll_layout.addWidget(chem_group)

        # --- Shared: labels group and operation parameters ---
        self._build_spec_labels_group(scroll_layout, labels_attr="spec_labels")
        self._build_spec_operation_group(scroll_layout)

        scroll_area.setWidget(scroll_widget)
        self.speciation_params_layout.addWidget(scroll_area)

        # Run / Stop buttons (mode-specific attribute names for backward compatibility)
        run_btn = QPushButton("▶ Generate Speciation Diagram")
        run_btn.setToolTip("Run the Speciation Diagram Generation")
        run_btn.setStyleSheet(ButtonStylesheets.RUN_BUTTON)
        run_btn.clicked.connect(self.run_compute_speciation_diagram)

        stop_button = QPushButton("⏹ Stop")
        stop_button.setToolTip("Stop the current operation.")
        stop_button.setStyleSheet(ButtonStylesheets.STOP_BUTTON)
        stop_button.setEnabled(False)
        stop_button.clicked.connect(self.stop_function)

        if pom_type == "IPA":
            self.ipa_comp_spec_run_btn = run_btn
            self.ipa_comp_spec_stop_btn = stop_button
        else:
            self.hpa_comp_spec_run_btn = run_btn
            self.hpa_comp_spec_stop_btn = stop_button

        self.speciation_params_layout.addWidget(run_btn)
        self.speciation_params_layout.addWidget(stop_button)

    def run_compute_speciation_diagram(self):
        # Start the speciation diagram generation process
        # ...
        self.clear_console()

        try:
            if self.pom_type == "IPA":
                from utilities.speciation_gui import speciation_ipa
                speciation_function = speciation_ipa
            elif self.pom_type == "HPA":
                from utilities.speciation_gui import speciation_hpa
                speciation_function = speciation_hpa
            else:
                QMessageBox.critical(self, "Error", "Invalid speciation type selected.")
                self.set_tab_error('speciation', 'spec_diag')
                return
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import speciation function: {str(e)}.")
            self.set_tab_error('speciation', 'spec_diag')
            return
        try:
            config_dict = self.get_gui_params()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to get GUI parameters: {str(e)}.")
            self.set_tab_error('speciation', 'spec_diag')
            return

        if not config_dict['Preparation']['output_path'] or not config_dict['Preparation']['POM_system']:
            QMessageBox.critical(self, "Error", "Output path or POM system not selected.")
            self.set_tab_error('speciation', 'spec_diag')
            return
        try:
            self.update_console("Started compute speciation", "speciation")
            self.set_tab_running('speciation', 'spec_diag')

            # Create and start function runner
            self.speciation_runner = POMSim_func_runner(speciation_function, config_dict)
            self.speciation_runner.finished.connect(self.on_function_finished)
            self.speciation_runner.error.connect(self.on_function_error)
            self.speciation_runner.stopped.connect(self.on_function_stopped)
            # Connect console output to update console
            self.speciation_runner.console_output.connect(lambda text: self.update_console(text, "speciation"))
            if self.pom_type == "IPA":
                self.ipa_comp_spec_run_btn.setEnabled(False)
                self.ipa_comp_spec_stop_btn.setEnabled(True)
            elif self.pom_type == "HPA":
                self.hpa_comp_spec_run_btn.setEnabled(False)
                self.hpa_comp_spec_stop_btn.setEnabled(True)
            self.speciation_runner.start()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to start speciation function: {str(e)}.")
            self.set_tab_error('speciation', 'spec_diag')

            if self.pom_type == "IPA":
                self.ipa_comp_spec_run_btn.setEnabled(True)
                self.ipa_comp_spec_stop_btn.setEnabled(False)
            elif self.pom_type == "HPA":
                self.hpa_comp_spec_run_btn.setEnabled(True)
                self.hpa_comp_spec_stop_btn.setEnabled(False)

    def create_phase_diagram_subtab(self, layout):
        """
        Create the Phase Diagram subtab with parameters for generating speciation diagrams.

        Args:
            layout (QVBoxLayout): The parent layout to add components to
        """
        # Speciation type selection

        # Parameters container
        self.phase_params_container = QWidget()
        self.phase_params_layout = QVBoxLayout(self.phase_params_container)
        layout.addWidget(self.phase_params_container)

        # Initialize parameters for Speciation
        self.update_phase_parameters()

        # Add console at bottom

    def _build_phase_io_group(self, layout):
        """Build the shared Input/Output File Management group for the phase diagram tab.

        Creates System, Output Path, Phase Diagram Directory Name, and Model Subset File
        fields. Registers System and Output Path as shared UI elements.
        Adds the resulting QGroupBox to *layout*.
        """
        system_group = QGroupBox("Input/Output File Management")
        system_layout = QVBoxLayout()

        system_row = QHBoxLayout()
        system_row.addWidget(QLabel("System:"))
        self.POM_system = QLineEdit()
        system_row.addWidget(self.POM_system)
        system_layout.addLayout(system_row)
        self.register_ui_element(self.POM_system, "POM_system")

        output_row = QHBoxLayout()
        output_row.addWidget(QLabel("Output Path:"))
        self.output_path = QLineEdit()
        output_row.addWidget(self.output_path)
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(lambda: self.browse_directory(self.output_path))
        output_row.addWidget(browse_btn)
        system_layout.addLayout(output_row)
        self.register_ui_element(self.output_path, "output_path")

        phase_dir_row = QHBoxLayout()
        phase_dir_row.addWidget(QLabel("Phase Diagram Directory Name:"))
        self.phase_dir = QLineEdit()
        phase_dir_row.addWidget(self.phase_dir)
        browse_btn_3 = QPushButton("Browse")
        browse_btn_3.clicked.connect(lambda: self.browse_directory(self.phase_dir))
        phase_dir_row.addWidget(browse_btn_3)
        system_layout.addLayout(phase_dir_row)
        self.register_ui_element(self.phase_dir, "phase_dir")

        model_subset_file_row = QHBoxLayout()
        model_subset_file_row.addWidget(QLabel("Model Subset File:"))
        self.model_subset_file = QLineEdit()
        model_subset_file_row.addWidget(self.model_subset_file)
        browse_btn2 = QPushButton("Browse")
        browse_btn2.clicked.connect(lambda: self.browse_file(self.model_subset_file))
        model_subset_file_row.addWidget(browse_btn2)
        system_layout.addLayout(model_subset_file_row)

        system_group.setLayout(system_layout)
        layout.addWidget(system_group)

    def _build_phase_ph_group(self, ph_attr_prefix="phase"):
        """Build a pH Range QGroupBox with min/max/step spinboxes for phase diagrams.

        Stores the created widgets as ``self.<ph_attr_prefix>_min_ph``,
        ``self.<ph_attr_prefix>_max_ph``, and ``self.<ph_attr_prefix>_step_ph``.

        Returns the configured QGroupBox.
        """
        ph_group = QGroupBox("pH Range")
        ph_layout = QVBoxLayout()

        min_ph_row = QHBoxLayout()
        min_ph_row.addWidget(QLabel("Min pH:"))
        min_ph = CustomDoubleSpinBox()
        min_ph.setRange(-10, 100)
        min_ph.setValue(0)
        min_ph_row.addWidget(min_ph)
        ph_layout.addLayout(min_ph_row)
        setattr(self, f"{ph_attr_prefix}_min_ph", min_ph)

        max_ph_row = QHBoxLayout()
        max_ph_row.addWidget(QLabel("Max pH:"))
        max_ph = CustomDoubleSpinBox()
        max_ph.setRange(-10, 100)
        max_ph.setValue(14)
        max_ph_row.addWidget(max_ph)
        ph_layout.addLayout(max_ph_row)
        setattr(self, f"{ph_attr_prefix}_max_ph", max_ph)

        step_ph_row = QHBoxLayout()
        step_ph_row.addWidget(QLabel("step pH:"))
        step_ph = CustomDoubleSpinBox()
        step_ph.setRange(0.01, 1)
        step_ph.setValue(0.5)
        step_ph.setDecimals(2)
        step_ph.setSingleStep(0.1)
        step_ph_row.addWidget(step_ph)
        ph_layout.addLayout(step_ph_row)
        setattr(self, f"{ph_attr_prefix}_step_ph", step_ph)

        ph_group.setLayout(ph_layout)
        return ph_group

    def _build_phase_labels_group(self, layout, labels_attr="phase_labels"):
        """Build the Species to Plot group with Labels File and Selected Labels fields for phase diagrams.

        Stores the labels file widget as ``self.labels_file`` and the selected labels
        widget as ``self.<labels_attr>``. Registers both as shared UI elements.
        Adds the resulting QGroupBox to *layout*.
        """
        plot_list_group = QGroupBox("Species to Plot")
        plot_list_layout = QVBoxLayout()
        plot_list_group.setLayout(plot_list_layout)

        labels_file_row = QHBoxLayout()
        labels_file_row.addWidget(QLabel("Labels File:"))
        self.labels_file = QLineEdit()
        labels_file_row.addWidget(self.labels_file)
        labels_file_browse_btn = QPushButton("Browse")
        labels_file_browse_btn.clicked.connect(lambda: self.browse_file(self.labels_file))
        self.register_ui_element(self.labels_file, "labels_file")
        labels_file_row.addWidget(labels_file_browse_btn)
        plot_list_layout.addLayout(labels_file_row)

        selected_labels_row = QHBoxLayout()
        selected_labels_row.addWidget(QLabel("Selected Labels:"))
        selected_labels = QLineEdit()
        selected_labels.setText("all")
        selected_labels.setReadOnly(False)
        selected_labels_row.addWidget(selected_labels)
        setattr(self, labels_attr, selected_labels)
        select_labels_btn = QPushButton("Select Labels")
        select_labels_btn.clicked.connect(
            lambda: self.load_and_select_labels(self.labels_file, getattr(self, labels_attr))
        )
        self.register_ui_element(selected_labels, labels_attr)
        selected_labels_row.addWidget(select_labels_btn)
        plot_list_layout.addLayout(selected_labels_row)

        layout.addWidget(plot_list_group)

    def _build_phase_operation_group(self, layout):
        """Build the Operation Parameters group with Cores and Batch Size spinboxes for phase diagrams.

        Stores widgets as ``self.phase_cores`` and ``self.phase_batch_size`` and
        registers them as shared UI elements. Adds the QGroupBox to *layout*.
        """
        params_group = QGroupBox("Operation Parameters")
        params_layout = QGridLayout()

        params_layout.addWidget(QLabel("Cores:"), 1, 0)
        self.phase_cores = QSpinBox()
        self.phase_cores.setRange(1, cpu_count() - 1)
        self.phase_cores.setValue(4)
        params_layout.addWidget(self.phase_cores, 1, 1)
        self.register_ui_element(self.phase_cores, "phase_cores")

        params_layout.addWidget(QLabel("Batch Size:"), 2, 0)
        self.phase_batch_size = QSpinBox()
        self.phase_batch_size.setRange(1, 1000)
        self.phase_batch_size.setValue(100)
        params_layout.addWidget(self.phase_batch_size, 2, 1)
        self.register_ui_element(self.phase_batch_size, "phase_batch_size")

        params_group.setLayout(params_layout)
        layout.addWidget(params_group)

    def update_phase_parameters(self):
        """
        Update the phase parameters UI based on the selected simulation type.
        """
        # Clear all existing widgets from the layout
        while self.phase_params_layout.count():
            child = self.phase_params_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

        pom_type = getattr(self, "global_mode", "IPA")

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)

        # --- Shared: I/O group ---
        self._build_phase_io_group(scroll_layout)

        # --- Chemical Parameters (mode-specific concentration + shared pH range) ---
        chem_group = QGroupBox("Chemical Parameters")
        chem_layout = QVBoxLayout()

        if pom_type == "IPA":
            min_c_row = QHBoxLayout()
            min_c_row.addWidget(QLabel("Min Concentration (mol/L):"))
            self.phase_min_c = CustomDoubleSpinBox()
            self.phase_min_c.setDecimals(6)
            self.phase_min_c.setRange(0.000001, 10)
            self.phase_min_c.setValue(0.001)
            self.phase_min_c.setSingleStep(0.1)
            min_c_row.addWidget(self.phase_min_c)
            chem_layout.addLayout(min_c_row)

            max_c_row = QHBoxLayout()
            max_c_row.addWidget(QLabel("Max Concentration (mol/L):"))
            self.phase_max_c = CustomDoubleSpinBox()
            self.phase_max_c.setDecimals(6)
            self.phase_max_c.setRange(0.000001, 10)
            self.phase_max_c.setValue(1.0)
            self.phase_max_c.setSingleStep(0.1)
            max_c_row.addWidget(self.phase_max_c)
            chem_layout.addLayout(max_c_row)

            num_c_row = QHBoxLayout()
            num_c_row.addWidget(QLabel("Number of Concentration points:"))
            self.phase_num_c = QSpinBox()
            self.phase_num_c.setMinimum(2)
            self.phase_num_c.setMaximum(1000)
            self.phase_num_c.setValue(10)
            num_c_row.addWidget(self.phase_num_c)
            chem_layout.addLayout(num_c_row)

        else:  # HPA
            min_ratio_row = QHBoxLayout()
            min_ratio_row.addWidget(QLabel("Min Metal/Heteroatom Ratio:"))
            self.phase_min_ratio = CustomDoubleSpinBox()
            self.phase_min_ratio.setDecimals(6)
            self.phase_min_ratio.setRange(0.000001, 500)
            self.phase_min_ratio.setValue(1.0)
            self.phase_min_ratio.setSingleStep(0.5)
            min_ratio_row.addWidget(self.phase_min_ratio)
            chem_layout.addLayout(min_ratio_row)

            max_ratio_row = QHBoxLayout()
            max_ratio_row.addWidget(QLabel("Max Metal/Heteroatom Ratio:"))
            self.phase_max_ratio = CustomDoubleSpinBox()
            self.phase_max_ratio.setDecimals(6)
            self.phase_max_ratio.setRange(0.000001, 500)
            self.phase_max_ratio.setValue(5.0)
            self.phase_max_ratio.setSingleStep(0.5)
            max_ratio_row.addWidget(self.phase_max_ratio)
            chem_layout.addLayout(max_ratio_row)

            num_ratio_row = QHBoxLayout()
            num_ratio_row.addWidget(QLabel("Number of Ratio points:"))
            self.phase_num_ratio = QSpinBox()
            self.phase_num_ratio.setMinimum(2)
            self.phase_num_ratio.setMaximum(1000)
            self.phase_num_ratio.setValue(10)
            num_ratio_row.addWidget(self.phase_num_ratio)
            chem_layout.addLayout(num_ratio_row)

            CX_row = QHBoxLayout()
            CX_row.addWidget(QLabel("Initial Heteroatom Concentration (mol/L):"))
            self.phase_CX = CustomDoubleSpinBox()
            self.phase_CX.setDecimals(6)
            self.phase_CX.setRange(0.000001, 10)
            self.phase_CX.setValue(0.1)
            self.phase_CX.setSingleStep(0.1)
            CX_row.addWidget(self.phase_CX)
            chem_layout.addLayout(CX_row)

        chem_layout.addWidget(self._build_spec_ph_group("phase"))
        chem_group.setLayout(chem_layout)
        scroll_layout.addWidget(chem_group)

        # --- Shared: labels group and operation parameters ---
        self._build_spec_labels_group(scroll_layout, labels_attr="spec_labels")
        self._build_spec_operation_group(scroll_layout)

        scroll_area.setWidget(scroll_widget)
        self.phase_params_layout.addWidget(scroll_area)

        # Run / Stop buttons (mode-specific attribute names for backward compatibility)
        run_btn = QPushButton("▶ Generate Phase Diagram")
        run_btn.setToolTip("Run the Phase Diagram Generation")
        run_btn.setStyleSheet(ButtonStylesheets.RUN_BUTTON)
        run_btn.clicked.connect(self.run_compute_phase_diagram)

        stop_button = QPushButton("⏹ Stop")
        stop_button.setToolTip("Stop the current operation.")
        stop_button.setStyleSheet(ButtonStylesheets.STOP_BUTTON)
        stop_button.setEnabled(False)
        stop_button.clicked.connect(self.stop_function)

        if pom_type == "IPA":
            self.ipa_phase_run_btn = run_btn
            self.ipa_phase_stop_btn = stop_button
        else:
            self.hpa_phase_run_btn = run_btn
            self.hpa_phase_stop_btn = stop_button

        self.phase_params_layout.addWidget(run_btn)
        self.phase_params_layout.addWidget(stop_button)

    def run_compute_phase_diagram(self):
        # Start the speciation diagram generation process
        # ...
        self.clear_console()

        try:
            if self.pom_type == "IPA":
                from utilities.speciation_gui import phase_ipa
                phase_function = phase_ipa
            elif self.pom_type == "HPA":
                from utilities.speciation_gui import phase_hpa
                phase_function = phase_hpa
            else:
                QMessageBox.critical(self, "Error", "Invalid phase type selected.")
                self.set_tab_error('speciation', 'phase_diag')
                return
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import speciation function: {str(e)}.")
            self.set_tab_error('speciation', 'phase_diag')
            return
        try:
            config_dict = self.get_gui_params()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to get GUI parameters: {str(e)}.")
            self.set_tab_error('speciation', 'phase_diag')
            return

        if not config_dict['Preparation']['output_path'] or not config_dict['Preparation']['POM_system']:
            QMessageBox.critical(self, "Error", "Output path or POM system not selected.")
            self.set_tab_error('speciation', 'phase_diag')
            return
        try:
            self.update_console("Started compute speciation", "speciation")
            self.set_tab_running('speciation', 'phase_diag')

            # Create and start function runner
            self.phase_runner = POMSim_func_runner(phase_function, config_dict)
            self.phase_runner.finished.connect(self.on_function_finished)
            self.phase_runner.error.connect(self.on_function_error)
            self.phase_runner.stopped.connect(self.on_function_stopped)
            # Connect console output to update console
            self.phase_runner.console_output.connect(lambda text: self.update_console(text, "speciation"))
            if self.pom_type == "IPA":
                self.ipa_phase_run_btn.setEnabled(False)
                self.ipa_phase_stop_btn.setEnabled(True)
            elif self.pom_type == "HPA":
                self.hpa_phase_run_btn.setEnabled(False)
                self.hpa_phase_stop_btn.setEnabled(True)
            self.phase_runner.start()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to start phase function: {str(e)}.")
            self.set_tab_error('speciation', 'phase_diag')

            if self.pom_type == "IPA":
                self.ipa_phase_run_btn.setEnabled(True)
                self.ipa_phase_stop_btn.setEnabled(False)
            elif self.pom_type == "HPA":
                self.hpa_phase_run_btn.setEnabled(True)
                self.hpa_phase_stop_btn.setEnabled(False)

    def create_clustering_tab(self):
        """
        Create the Clustering tab with three subtabs for different clustering analyses.

        This method initializes the Clustering tab in the GUI, which contains three subtabs:
        1. SM_Clusterization: For performing basic clustering of simulation models
        2. Cluster_Model_Selection: For selecting optimal clustering models based on metrics
        3. Clust_Boxplot_Filtering: For filtering clusters using boxplot statistics

        The tab includes a tabbed interface for the three clustering types and a console output
        area at the bottom for displaying execution logs related to clustering operations.

        Parameters:
            self: The POMSimulatorGUI instance containing the tab widget where
              this tab will be added.

        Returns:
            None: This method doesn't return a value but modifies the GUI by adding
              a new tab to the main tab widget.
        """
        tab = QWidget()
        self.tabs.addTab(tab, "Clustering")
        layout = QVBoxLayout(tab)

        # Create subtabs
        subtabs = QTabWidget()
        layout.addWidget(subtabs)

        # Create separate color dict selector instances for each subtab
        self.clustering_col_dict = ColorDictSelector()
        self.register_ui_element(self.clustering_col_dict, "col_dict")

        self.filtering_col_dict = ColorDictSelector()
        self.register_ui_element(self.filtering_col_dict, "col_dict")

        # Connect the two selectors so they stay in sync
        # When clustering_col_dict changes, update filtering_col_dict
        self.clustering_col_dict.col_dict_changed.connect(
            self.filtering_col_dict.col_dict_combo.setCurrentText
        )
        # When filtering_col_dict changes, update clustering_col_dict
        self.filtering_col_dict.col_dict_changed.connect(
            self.clustering_col_dict.col_dict_combo.setCurrentText
        )

        # clustering subtab
        self.clustering_widget = QWidget()
        self.clustering_layout = QVBoxLayout(self.clustering_widget)
        self.create_clustering_subtab(self.clustering_layout)
        subtabs.addTab(self.clustering_widget, "Clustering")

        # selection subtab
        self.selection_widget = QWidget()
        self.selection_layout = QVBoxLayout(self.selection_widget)
        self.create_selection_subtab(self.selection_layout)
        subtabs.addTab(self.selection_widget, "Model Selection")

        # filtering subtab
        self.filtering_widget = QWidget()
        self.filtering_layout = QVBoxLayout(self.filtering_widget)
        self.create_filtering_subtab(self.filtering_layout)
        subtabs.addTab(self.filtering_widget, "Boxplot Filtering")

    def load_labels_from_file(self, file_path):
        """
        Load species labels from a labels file.

        Args:
            file_path (str): Path to the labels file

        Returns:
            dict: Dictionary with species names as keys, or empty dict if file not found
        """
        try:
            labels_dict = {}
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line:
                        labels_dict[line] = False  # Initialize as unchecked

            # print(f"Loaded {len(labels_dict)} species from {file_path}")
            return labels_dict

        except Exception as e:
            # print(f"Error loading labels file: {e}")
            QMessageBox.warning(self, "Error Loading Labels", f"Failed to load labels file: {str(e)}")
            return {}

    def load_and_select_labels(self, labels_file, selected_labels_widget):
        """
        Load labels from file and open selection dialog.

        Args:
            labels_file (QLineEdit): The QLineEdit containing the labels file path
            selected_labels_widget (QLineEdit): The QLineEdit to store selected labels
        """

        if not labels_file.text():
            QMessageBox.warning(self, "Warning", "Please select a labels file first.")
            return

        labels_file_path = script_dir + labels_file.text().split("..")[1]
        file_path = os.path.abspath(labels_file_path)

        if not file_path:
            QMessageBox.warning(self, "Warning", "Please select a labels file first.")
            return

        if not os.path.exists(file_path):
            QMessageBox.critical(self, "Error", f"Labels file not found: {file_path}")
            return

        try:
            # Load original labels from file
            original_labels_list = list(self.load_labels_from_file(file_path).keys())

            if not original_labels_list:
                QMessageBox.warning(self, "Warning", "No labels found in the selected file.")
                return

            # Open selection dialog
            dialog = LabelSelectionDialog(original_labels_list, self)
            if dialog.exec_() == QDialog.Accepted:
                selected_labels = dialog.get_selected_labels()

                # Compare selected labels with original list
                # If they are exactly the same, use empty string (means "all")
                if sorted(selected_labels) == sorted(original_labels_list):
                    labels_to_store = "all"
                else:
                    labels_to_store = ",".join(selected_labels)

                # Store the selected labels in the widget
                selected_labels_widget.setText(labels_to_store)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load labels: {str(e)}")

    def _build_clustering_io_group(self, layout, include_features_file=True, include_m_idx=False):
        """Build the shared Input Output File Managing group for clustering subtabs.

        Creates System, Output path, Cluster directory, and NPZ File fields.
        Optionally adds a Features File row (for clustering/selection subtabs) or
        an m_idx ComboBox (for the filtering subtab).

        Registers all fields as shared UI elements and adds the QGroupBox to *layout*.
        """
        system_group = QGroupBox("Input Output File Managing")
        system_layout = QVBoxLayout()

        system_row = QHBoxLayout()
        system_row.addWidget(QLabel("System:"))
        self.POM_system = QLineEdit()
        system_row.addWidget(self.POM_system)
        system_layout.addLayout(system_row)
        self.register_ui_element(self.POM_system, "POM_system")

        output_row = QHBoxLayout()
        output_row.addWidget(QLabel("Output path:"))
        self.output_path = QLineEdit()
        output_row.addWidget(self.output_path)
        browse_output_btn = QPushButton("Browse")
        browse_output_btn.clicked.connect(lambda: self.browse_directory(self.output_path))
        output_row.addWidget(browse_output_btn)
        system_layout.addLayout(output_row)
        self.register_ui_element(self.output_path, "output_path")

        cluster_dir_row = QHBoxLayout()
        cluster_dir_row.addWidget(QLabel("Cluster directory:"))
        self.cluster_dir = QLineEdit()
        cluster_dir_row.addWidget(self.cluster_dir)
        dir_browse_btn = QPushButton("Browse")
        dir_browse_btn.clicked.connect(lambda: self.browse_directory(self.cluster_dir))
        cluster_dir_row.addWidget(dir_browse_btn)
        system_layout.addLayout(cluster_dir_row)
        self.register_ui_element(self.cluster_dir, "cluster_dir")

        npz_cluster_file_row = QHBoxLayout()
        npz_cluster_file_row.addWidget(QLabel("NPZ File:"))
        self.npz_cluster_file = QLineEdit()
        npz_cluster_file_row.addWidget(self.npz_cluster_file)
        npz_browse_btn = QPushButton("Browse")
        npz_browse_btn.clicked.connect(lambda: self.browse_file(self.npz_cluster_file, "NPZ Files (*.npz)"))
        npz_cluster_file_row.addWidget(npz_browse_btn)
        system_layout.addLayout(npz_cluster_file_row)
        self.register_ui_element(self.npz_cluster_file, "npz_cluster_file")

        if include_features_file:
            features_file_row = QHBoxLayout()
            features_file_row.addWidget(QLabel("Features File:"))
            self.features_file = QLineEdit()
            features_file_row.addWidget(self.features_file)
            features_browse_btn = QPushButton("Browse")
            features_browse_btn.clicked.connect(lambda: self.browse_file(self.features_file, "CSV Files (*.csv)"))
            features_file_row.addWidget(features_browse_btn)
            system_layout.addLayout(features_file_row)
            self.register_ui_element(self.features_file, "features_file")

        if include_m_idx:
            m_idx_row = QHBoxLayout()
            m_idx_row.addWidget(QLabel("m_idx:"))
            self.m_idx = QComboBox()
            self.m_idx.addItems(["0", "1"])
            m_idx_row.addWidget(self.m_idx)
            system_layout.addLayout(m_idx_row)

        system_group.setLayout(system_layout)
        layout.addWidget(system_group)
        return system_group, system_layout

    def create_clustering_subtab(self, layout):

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)

        # Shared file management group (System, Output, Cluster dir, NPZ, Features)
        self._build_clustering_io_group(scroll_layout, include_features_file=True, include_m_idx=False)

        # Operation parameters
        operation_group = QGroupBox("Operation Parameters")
        operation_layout = QVBoxLayout()

        n_clusters_row = QHBoxLayout()
        n_clusters_row.addWidget(QLabel("Number of Clusters:"))
        self.n_clusters = CustomDoubleSpinBox()
        self.n_clusters.setMinimum(2)
        self.n_clusters.valueChanged.connect(self.update_cluster_selection)
        n_clusters_row.addWidget(self.n_clusters)
        operation_layout.addLayout(n_clusters_row)

        normalize_feats_row = QHBoxLayout()
        normalize_feats_row.addWidget(QLabel("Normalize Features:"))
        self.normalize_feats = QComboBox()
        self.normalize_feats.addItems(["True", "False"])
        normalize_feats_row.addWidget(self.normalize_feats)
        operation_layout.addLayout(normalize_feats_row)

        features_group = QGroupBox("Feature Selection")
        feature_selection_layout = QGridLayout()
        self.sel_features = {}
        features_list = clustering_features

        ncols = len(features_list)
        for jj, feature in enumerate(features_list):
            row = jj // ncols
            col = jj % ncols

            checkbox = QCheckBox(feature)
            checkbox.setChecked(True)
            self.sel_features[feature] = checkbox
            feature_selection_layout.addWidget(checkbox, row, col)

        features_group.setLayout(feature_selection_layout)
        operation_layout.addWidget(features_group)

        operation_group.setLayout(operation_layout)
        scroll_layout.addWidget(operation_group)

        # Add plot list (species selection) group
        plot_list_group = QGroupBox("Species to Plot")
        plot_list_layout = QVBoxLayout()
        plot_list_group.setLayout(plot_list_layout)

        # Labels file row
        labels_file_row = QHBoxLayout()
        labels_file_row.addWidget(QLabel("Labels File:"))
        self.labels_file = QLineEdit()
        labels_file_row.addWidget(self.labels_file)
        labels_file_browse_btn = QPushButton("Browse")
        labels_file_browse_btn.clicked.connect(
            lambda: self.browse_file(self.labels_file)
        )
        self.register_ui_element(self.labels_file, "labels_file")

        labels_file_row.addWidget(labels_file_browse_btn)
        plot_list_layout.addLayout(labels_file_row)

        # Selected labels row
        selected_labels_row = QHBoxLayout()
        selected_labels_row.addWidget(QLabel("Selected Labels:"))
        self.plot_list = QLineEdit()
        self.plot_list.setText("all")
        self.plot_list.setReadOnly(False)
        selected_labels_row.addWidget(self.plot_list)
        select_labels_btn = QPushButton("Select Labels")
        select_labels_btn.clicked.connect(
            lambda: self.load_and_select_labels(self.labels_file, self.plot_list)
        )
        self.register_ui_element(self.plot_list, "plot_list")
        selected_labels_row.addWidget(select_labels_btn)
        plot_list_layout.addLayout(selected_labels_row)

        scroll_layout.addWidget(plot_list_group)

        # Add color dictionary selector
        color_dict_group = QGroupBox("Color Dictionary Settings")
        color_dict_layout = QVBoxLayout()
        color_dict_layout.addWidget(self.clustering_col_dict)
        color_dict_group.setLayout(color_dict_layout)
        scroll_layout.addWidget(color_dict_group)

        # Add stretch to push everything to the top
        scroll_layout.addStretch()

        scroll_area.setWidget(scroll_widget)
        layout.addWidget(scroll_area)

        # Add buttons below scroll area
        run_btn = QPushButton("▶ Run Clustering")
        run_btn.setToolTip("Run clustering operation")
        run_btn.setStyleSheet(ButtonStylesheets.RUN_BUTTON)
        run_btn.clicked.connect(self.run_compute_clustering)
        self.cluster_run_btn = run_btn
        layout.addWidget(self.cluster_run_btn)

        stop_btn = QPushButton("⏹ Stop")
        stop_btn.setToolTip("Stop the current function.")
        stop_btn.setStyleSheet(ButtonStylesheets.STOP_BUTTON)
        stop_btn.clicked.connect(self.stop_function)
        stop_btn.setEnabled(False)
        self.cluster_stop_btn = stop_btn
        layout.addWidget(self.cluster_stop_btn)

        return

    def run_compute_clustering(self):
        # Start the speciation diagram generation process
        # ...
        self.clear_console()

        try:
            from utilities.clustering_gui import SM_clustering

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import speciation function: {str(e)}.")
            self.set_tab_error('clustering', 'clustering')
            return
        try:
            config_dict = self.get_gui_params()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to get GUI parameters: {str(e)}.")
            self.set_tab_error('clustering', 'clustering')
            return

        if not config_dict['Preparation']['output_path'] or not config_dict['Preparation']['POM_system']:
            QMessageBox.critical(self, "Error", "Output path or POM system not selected.")
            self.set_tab_error('clustering', 'clustering')
            return
        try:
            self.update_console("Started compute speciation", "speciation")
            self.set_tab_running('clustering', 'clustering')

            # Create and start function runner
            self.cluster_runner = POMSim_func_runner(SM_clustering, config_dict)
            self.cluster_runner.finished.connect(self.on_function_finished)
            self.cluster_runner.error.connect(self.on_function_error)
            self.cluster_runner.stopped.connect(self.on_function_stopped)
            # Connect console output to update console
            self.cluster_runner.console_output.connect(lambda text: self.update_console(text, "clustering"))

            self.cluster_run_btn.setEnabled(False)
            self.cluster_stop_btn.setEnabled(True)

            self.cluster_runner.start()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to start phase function: {str(e)}.")
            self.set_tab_error('clustering', 'clustering')

            self.cluster_run_btn.setEnabled(True)
            self.cluster_stop_btn.setEnabled(False)

    def update_cluster_selection(self):
        """
        Updates the cluster selection checkboxes in the Model Selection tab
        when the number of clusters is changed in the Clustering tab.
        """
        # Get the current number of clusters
        n_clusters = int(self.n_clusters.value())

        # Find the layout in the selection tab
        for i in range(self.selection_widget.layout().count()):
            item = self.selection_widget.layout().itemAt(i)
            if isinstance(item.widget(), QScrollArea):
                scroll_area = item.widget()
                scroll_content = scroll_area.widget()

                # Find the system group
                for j in range(scroll_content.layout().count()):
                    group_item = scroll_content.layout().itemAt(j)
                    if isinstance(group_item.widget(),
                                  QGroupBox) and group_item.widget().title() == "Input Output File Managing":
                        system_group = group_item.widget()

                        # Find the Feature Selection group within the system group
                        for k in range(system_group.layout().count()):
                            child_item = system_group.layout().itemAt(k)
                            if isinstance(child_item.widget(),
                                          QGroupBox) and child_item.widget().title() == "Cluster Selection":
                                sel_group_idx_group = child_item.widget()
                                sel_group_idx_layout = sel_group_idx_group.layout()

                                # Clear existing checkboxes
                                self.selected_group_idx = {}
                                while sel_group_idx_layout.count():
                                    item = sel_group_idx_layout.takeAt(0)
                                    if item.widget():
                                        item.widget().deleteLater()

                                # Add new checkboxes
                                ncols = 5
                                for kk in range(n_clusters):
                                    row = kk // ncols
                                    col = kk % ncols

                                    checkbox = QCheckBox(str(kk))
                                    checkbox.setChecked(False)
                                    self.selected_group_idx[str(kk)] = checkbox
                                    sel_group_idx_layout.addWidget(checkbox, row, col)

                                # Update the layout
                                sel_group_idx_group.updateGeometry()
                                break

        # Force update of the UI
        self.selection_widget.update()

    def create_selection_subtab(self, layout):

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)

        # Shared file management group (System, Output, Cluster dir, NPZ, Features)
        system_group, system_layout = self._build_clustering_io_group(scroll_layout, include_features_file=True,
                                                                      include_m_idx=False)

        sel_group_idx_group = QGroupBox("Cluster Selection")
        sel_group_idx_layout = QGridLayout()
        self.selected_group_idx = {}
        n_clusters = int(self.n_clusters.value())

        ncols = 5
        for kk in range(n_clusters):
            row = kk // ncols
            col = kk % ncols

            checkbox = QCheckBox(str(kk))
            checkbox.setChecked(False)
            self.selected_group_idx[str(kk)] = checkbox
            sel_group_idx_layout.addWidget(checkbox, row, col)

        sel_group_idx_group.setLayout(sel_group_idx_layout)
        system_layout.addWidget(sel_group_idx_group)

        system_group.setLayout(system_layout)
        scroll_layout.addWidget(system_group)

        scroll_widget.setLayout(scroll_layout)
        scroll_area.setWidget(scroll_widget)
        self.selection_layout.addWidget(scroll_area)

        run_btn = QPushButton("Run cluster selection")
        run_btn.setToolTip("Run clustering operation")
        run_btn.setStyleSheet(ButtonStylesheets.RUN_BUTTON)

        run_btn.clicked.connect(self.run_compute_selection)
        self.selection_run_btn = run_btn  # Store reference to enable/disable it
        self.selection_layout.addWidget(self.selection_run_btn)

        stop_btn = QPushButton("⏹ Stop")
        stop_btn.setToolTip("Stop the current simulation.")
        stop_btn.setStyleSheet(ButtonStylesheets.STOP_BUTTON)
        stop_btn.clicked.connect(self.stop_function)
        stop_btn.setEnabled(False)  # Initially disabled
        self.selection_stop_btn = stop_btn  # Store reference to enable/disable it
        self.selection_layout.addWidget(self.selection_stop_btn)

        return

    def run_compute_selection(self):

        self.clear_console()

        try:
            from utilities.clustering_gui import cluster_selection

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import speciation function: {str(e)}.")
            self.set_tab_error('clustering', 'selection')
            return
        try:
            config_dict = self.get_gui_params()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to get GUI parameters: {str(e)}.")
            self.set_tab_error('clustering', 'selection')
            return

        if not config_dict['Preparation']['output_path'] or not config_dict['Preparation']['POM_system']:
            QMessageBox.critical(self, "Error", "Output path or POM system not selected.")
            self.set_tab_error('clustering', 'selection')
            return
        try:
            self.update_console("Started compute speciation", "speciation")
            self.set_tab_running('clustering', 'selection')

            # Create and start function runner
            self.selection_runner = POMSim_func_runner(cluster_selection, config_dict)
            self.selection_runner.finished.connect(self.on_function_finished)
            self.selection_runner.error.connect(self.on_function_error)
            self.selection_runner.stopped.connect(self.on_function_stopped)
            # Connect console output to update console
            self.selection_runner.console_output.connect(lambda text: self.update_console(text, "clustering"))

            self.selection_run_btn.setEnabled(False)
            self.selection_stop_btn.setEnabled(True)

            self.selection_runner.start()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to start phase function: {str(e)}.")
            self.set_tab_error('clustering', 'selection')

            self.selection_run_btn.setEnabled(True)
            self.selection_stop_btn.setEnabled(False)

    def create_filtering_subtab(self, layout):

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)

        # Shared file management group (System, Output, Cluster dir, NPZ, m_idx)
        self._build_clustering_io_group(scroll_layout, include_features_file=False, include_m_idx=True)

        # Add plot list (species selection) group
        plot_list_group = QGroupBox("Species to Plot")
        plot_list_layout = QVBoxLayout()
        plot_list_group.setLayout(plot_list_layout)

        # Labels file row
        labels_file_row = QHBoxLayout()
        labels_file_row.addWidget(QLabel("Labels File:"))
        self.labels_file = QLineEdit()
        labels_file_row.addWidget(self.labels_file)
        labels_file_browse_btn = QPushButton("Browse")
        labels_file_browse_btn.clicked.connect(
            lambda: self.browse_file(self.labels_file)
        )
        self.register_ui_element(self.labels_file, "labels_file")

        labels_file_row.addWidget(labels_file_browse_btn)
        plot_list_layout.addLayout(labels_file_row)

        # Selected plotting labels row
        plot_list_labels_row = QHBoxLayout()
        plot_list_labels_row.addWidget(QLabel("Selected plotting labels:"))
        self.plot_list = QLineEdit()
        self.plot_list.setText("all")
        self.plot_list.setReadOnly(False)
        plot_list_labels_row.addWidget(self.plot_list)
        select_labels_btn = QPushButton("Select Labels")
        select_labels_btn.clicked.connect(
            lambda: self.load_and_select_labels(self.labels_file, self.plot_list)
        )
        self.register_ui_element(self.plot_list, "plot_list")
        plot_list_labels_row.addWidget(select_labels_btn)
        plot_list_layout.addLayout(plot_list_labels_row)

        # Selected boxplot labels row
        boxplot_list_labels_row = QHBoxLayout()
        boxplot_list_labels_row.addWidget(QLabel("Selected filtering labels:"))
        self.boxplot_list = QLineEdit()
        self.boxplot_list.setText("all")
        self.boxplot_list.setReadOnly(False)
        boxplot_list_labels_row.addWidget(self.boxplot_list)
        select_labels_btn = QPushButton("Select Labels")
        select_labels_btn.clicked.connect(
            lambda: self.load_and_select_labels(self.labels_file, self.boxplot_list)
        )
        boxplot_list_labels_row.addWidget(select_labels_btn)
        plot_list_layout.addLayout(boxplot_list_labels_row)

        scroll_layout.addWidget(plot_list_group)

        color_dict_group = QGroupBox("Color Dictionary Settings")
        color_dict_layout = QVBoxLayout()
        color_dict_layout.addWidget(self.filtering_col_dict)
        color_dict_group.setLayout(color_dict_layout)
        scroll_layout.addWidget(color_dict_group)

        # Add stretch to push everything to the top
        scroll_layout.addStretch()

        scroll_widget.setLayout(scroll_layout)
        scroll_area.setWidget(scroll_widget)
        self.filtering_layout.addWidget(scroll_area)

        run_btn = QPushButton("Run Filtering")
        run_btn.setToolTip("Run filtering operation")
        run_btn.setStyleSheet(ButtonStylesheets.RUN_BUTTON)

        run_btn.clicked.connect(self.run_compute_filtering)
        self.filtering_run_btn = run_btn  # Store reference to enable/disable it
        self.filtering_layout.addWidget(self.filtering_run_btn)

        stop_btn = QPushButton("⏹ Stop")
        stop_btn.setToolTip("Stop the current simulation.")
        stop_btn.setStyleSheet(ButtonStylesheets.STOP_BUTTON)
        stop_btn.clicked.connect(self.stop_function)
        stop_btn.setEnabled(False)  # Initially disabled
        self.filtering_stop_btn = stop_btn  # Store reference to enable/disable it
        self.filtering_layout.addWidget(self.filtering_stop_btn)

        return

    def run_compute_filtering(self):

        self.clear_console()

        try:
            from utilities.clustering_gui import clust_filtering

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import speciation function: {str(e)}.")
            self.set_tab_error('clustering', 'filtering')
            return
        try:
            config_dict = self.get_gui_params()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to get GUI parameters: {str(e)}.")
            self.set_tab_error('clustering', 'filtering')
            return

        if not config_dict['Preparation']['output_path'] or not config_dict['Preparation']['POM_system']:
            QMessageBox.critical(self, "Error", "Output path or POM system not selected.")
            self.set_tab_error('clustering', 'filtering')
            return
        try:
            self.update_console("Started filtering", "clustering")
            self.set_tab_running('clustering', 'filtering')

            # Create and start function runner
            self.filtering_runner = POMSim_func_runner(clust_filtering, config_dict)
            self.filtering_runner.finished.connect(self.on_function_finished)
            self.filtering_runner.error.connect(self.on_function_error)
            self.filtering_runner.stopped.connect(self.on_function_stopped)
            # Connect console output to update console
            self.filtering_runner.console_output.connect(lambda text: self.update_console(text, "clustering"))

            self.filtering_run_btn.setEnabled(False)
            self.filtering_stop_btn.setEnabled(True)

            self.filtering_runner.start()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to start phase function: {str(e)}.")
            self.set_tab_error('clustering', 'filtering')

            self.filtering_run_btn.setEnabled(True)
            self.filtering_stop_btn.setEnabled(False)

    def create_plotting_tab(self):
        """
        Create the Plotting tab with two subtabs for different plotting types.

        This method initializes the Plotting tab in the GUI, which contains two subtabs:
        1. Plot Speciation Diagram: For plotting speciation diagrams from NPZ files
        2. Plot Phase Diagram: For plotting phase diagrams from phase calculation results

        The tab includes a tabbed interface for the two plotting types and uses the
        standard threaded execution infrastructure for non-blocking operations.
        """
        tab = QWidget()
        self.tabs.addTab(tab, "Plotting")
        layout = QVBoxLayout(tab)

        # Create subtabs
        subtabs = QTabWidget()
        layout.addWidget(subtabs)

        # Plot Speciation Diagram subtab
        self.plot_spec_widget = QWidget()
        self.plot_spec_layout = QVBoxLayout(self.plot_spec_widget)
        self.create_plot_speciation_subtab(self.plot_spec_layout)
        subtabs.addTab(self.plot_spec_widget, "Plot Speciation Diagram")

        # Plot Phase Diagram subtab
        self.plot_phase_widget = QWidget()
        self.plot_phase_layout = QVBoxLayout(self.plot_phase_widget)
        self.create_plot_phase_subtab(self.plot_phase_layout)
        subtabs.addTab(self.plot_phase_widget, "Plot Phase Diagram")

    def create_plot_speciation_subtab(self, layout):
        """
        Create the Plot Speciation Diagram subtab with parameters for plotting speciation diagrams.

        Args:
            layout (QVBoxLayout): The parent layout to add components to
        """
        # Parameters container
        self.plot_spec_params_container = QWidget()
        self.plot_spec_params_layout = QVBoxLayout(self.plot_spec_params_container)
        layout.addWidget(self.plot_spec_params_container)

        # Initialize parameters for Plot Speciation
        self.update_plot_speciation_parameters()

    def update_plot_speciation_parameters(self):
        """
        Update the plot speciation parameters UI based on the selected simulation type.
        """
        # Clear all existing widgets from the layout
        while self.plot_spec_params_layout.count():
            child = self.plot_spec_params_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

        pom_type = getattr(self, "global_mode", "IPA")

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)

        # --- Input/Output File Management ---
        system_group = QGroupBox("Input/Output File Management")
        system_layout = QVBoxLayout()

        # System row
        system_row = QHBoxLayout()
        system_row.addWidget(QLabel("System:"))
        self.plot_spec_POM_system = QLineEdit()
        self.plot_spec_POM_system.setPlaceholderText("Enter POM system (e.g., As, W, PMo)")
        system_row.addWidget(self.plot_spec_POM_system)
        system_layout.addLayout(system_row)
        self.register_ui_element(self.plot_spec_POM_system, "POM_system")

        # Output path row
        output_row = QHBoxLayout()
        output_row.addWidget(QLabel("Output Path:"))
        self.plot_spec_output_path = QLineEdit()
        self.plot_spec_output_path.setPlaceholderText("Select output directory")
        output_row.addWidget(self.plot_spec_output_path)
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(lambda: self.browse_directory(self.plot_spec_output_path))
        output_row.addWidget(browse_btn)
        system_layout.addLayout(output_row)
        self.register_ui_element(self.plot_spec_output_path, "output_path")

        # NPZ file row
        npz_row = QHBoxLayout()
        npz_row.addWidget(QLabel("NPZ File:"))
        self.npz_file = self.get_or_create_widget('npz_file', QLineEdit)
        npz_row.addWidget(self.npz_file)
        npz_browse_btn = QPushButton("Browse")
        npz_browse_btn.clicked.connect(lambda: self.browse_file(self.npz_file, "NPZ Files (*.npz)"))
        npz_row.addWidget(npz_browse_btn)
        system_layout.addLayout(npz_row)
        self.register_ui_element(self.npz_file, "npz_file")

        # m_idx row (metal index for plotting)
        m_idx_row = QHBoxLayout()
        m_idx_row.addWidget(QLabel("Metal Index (m_idx):"))
        self.m_idx = self.get_or_create_widget('m_idx', QComboBox)
        if pom_type == "IPA":
            self.m_idx.clear()
            self.m_idx.addItems(["0"])
        else:  # HPA
            self.m_idx.clear()
            self.m_idx.addItems(["0", "1"])
        m_idx_row.addWidget(self.m_idx)
        system_layout.addLayout(m_idx_row)
        self.register_ui_element(self.m_idx, "m_idx")

        system_group.setLayout(system_layout)
        scroll_layout.addWidget(system_group)

        # --- Species to Plot ---
        plot_list_group = QGroupBox("Species to Plot")
        plot_list_layout = QVBoxLayout()
        plot_list_group.setLayout(plot_list_layout)

        # Labels file row
        labels_file_row = QHBoxLayout()
        labels_file_row.addWidget(QLabel("Labels File:"))
        # Create separate widget for speciation subtab
        if not hasattr(self, 'spec_labels_file'):
            self.spec_labels_file = QLineEdit()
            self.register_ui_element(self.spec_labels_file, "labels_file")
        labels_file_row.addWidget(self.spec_labels_file)
        labels_file_browse_btn = QPushButton("Browse")
        labels_file_browse_btn.clicked.connect(lambda: self.browse_file(self.spec_labels_file))
        labels_file_row.addWidget(labels_file_browse_btn)
        plot_list_layout.addLayout(labels_file_row)

        # Selected labels row
        selected_labels_row = QHBoxLayout()
        selected_labels_row.addWidget(QLabel("Selected Labels:"))
        # Create separate widget for speciation subtab
        if not hasattr(self, 'spec_plot_list'):
            self.spec_plot_list = QLineEdit()
            self.spec_plot_list.setText("all")
            self.spec_plot_list.setReadOnly(False)
            self.register_ui_element(self.spec_plot_list, "plot_list")
        selected_labels_row.addWidget(self.spec_plot_list)
        select_labels_btn = QPushButton("Select Labels")
        select_labels_btn.clicked.connect(
            lambda: self.load_and_select_labels(self.spec_labels_file, self.spec_plot_list)
        )
        selected_labels_row.addWidget(select_labels_btn)
        plot_list_layout.addLayout(selected_labels_row)

        scroll_layout.addWidget(plot_list_group)

        # --- Color Dictionary ---
        color_dict_group = QGroupBox("Color Dictionary")
        color_dict_layout = QVBoxLayout()
        
        # Create or get existing color dict selector
        if not hasattr(self, 'plot_spec_col_dict'):
            self.plot_spec_col_dict = ColorDictSelector()
            self.register_ui_element(self.plot_spec_col_dict, "col_dict")
        
        color_dict_layout.addWidget(self.plot_spec_col_dict)
        color_dict_group.setLayout(color_dict_layout)
        scroll_layout.addWidget(color_dict_group)

        scroll_area.setWidget(scroll_widget)
        self.plot_spec_params_layout.addWidget(scroll_area)

        # Run / Stop buttons
        run_btn = QPushButton("▶ Generate Speciation Plot")
        run_btn.setToolTip("Generate the Speciation Diagram Plot")
        run_btn.setStyleSheet(ButtonStylesheets.RUN_BUTTON)
        run_btn.clicked.connect(self.run_plot_speciation)

        stop_button = QPushButton("⏹ Stop")
        stop_button.setToolTip("Stop the current operation.")
        stop_button.setStyleSheet(ButtonStylesheets.STOP_BUTTON)
        stop_button.setEnabled(False)
        stop_button.clicked.connect(self.stop_function)

        self.plot_spec_run_btn = run_btn
        self.plot_spec_stop_btn = stop_button

        self.plot_spec_params_layout.addWidget(run_btn)
        self.plot_spec_params_layout.addWidget(stop_button)

    def create_plot_phase_subtab(self, layout):
        """
        Create the Plot Phase Diagram subtab with parameters for plotting phase diagrams.

        Args:
            layout (QVBoxLayout): The parent layout to add components to
        """
        # Parameters container
        self.plot_phase_params_container = QWidget()
        self.plot_phase_params_layout = QVBoxLayout(self.plot_phase_params_container)
        layout.addWidget(self.plot_phase_params_container)

        # Initialize parameters for Plot Phase
        self.update_plot_phase_parameters()

    def update_plot_phase_parameters(self):
        """
        Update the plot phase parameters UI based on the selected simulation type.
        """
        # Clear all existing widgets from the layout
        while self.plot_phase_params_layout.count():
            child = self.plot_phase_params_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

        pom_type = getattr(self, "global_mode", "IPA")

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout(scroll_widget)

        # --- Input/Output File Management ---
        system_group = QGroupBox("Input/Output File Management")
        system_layout = QVBoxLayout()

        # System row
        system_row = QHBoxLayout()
        system_row.addWidget(QLabel("System:"))
        self.plot_phase_POM_system = QLineEdit()
        self.plot_phase_POM_system.setPlaceholderText("Enter POM system (e.g., As, W, PMo)")
        system_row.addWidget(self.plot_phase_POM_system)
        system_layout.addLayout(system_row)
        self.register_ui_element(self.plot_phase_POM_system, "POM_system")

        # Output path row
        output_row = QHBoxLayout()
        output_row.addWidget(QLabel("Output Path:"))
        self.plot_phase_output_path = QLineEdit()
        self.plot_phase_output_path.setPlaceholderText("Select output directory")
        output_row.addWidget(self.plot_phase_output_path)
        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(lambda: self.browse_directory(self.plot_phase_output_path))
        output_row.addWidget(browse_btn)
        system_layout.addLayout(output_row)
        self.register_ui_element(self.plot_phase_output_path, "output_path")

        # Phase diagram directory row
        plot_phase_dir_row = QHBoxLayout()
        plot_phase_dir_row.addWidget(QLabel("Phase Diagram Directory Name:"))
        self.plot_phase_dir = QLineEdit()
        plot_phase_dir_row.addWidget(self.plot_phase_dir)
        system_layout.addLayout(plot_phase_dir_row)
        self.register_ui_element(self.plot_phase_dir, "phase_dir")

        system_group.setLayout(system_layout)
        scroll_layout.addWidget(system_group)

        # --- Species to Plot ---
        plot_list_group = QGroupBox("Species to Plot")
        plot_list_layout = QVBoxLayout()
        plot_list_group.setLayout(plot_list_layout)

        # Labels file row
        labels_file_row = QHBoxLayout()
        labels_file_row.addWidget(QLabel("Labels File:"))
        # Create separate widget for phase subtab
        if not hasattr(self, 'phase_labels_file'):
            self.phase_labels_file = QLineEdit()
            self.register_ui_element(self.phase_labels_file, "labels_file")
        labels_file_row.addWidget(self.phase_labels_file)
        labels_file_browse_btn = QPushButton("Browse")
        labels_file_browse_btn.clicked.connect(lambda: self.browse_file(self.phase_labels_file))
        labels_file_row.addWidget(labels_file_browse_btn)
        plot_list_layout.addLayout(labels_file_row)

        # Selected labels row
        selected_labels_row = QHBoxLayout()
        selected_labels_row.addWidget(QLabel("Selected Labels:"))
        # Create separate widget for phase subtab
        if not hasattr(self, 'phase_plot_list'):
            self.phase_plot_list = QLineEdit()
            self.phase_plot_list.setText("all")
            self.phase_plot_list.setReadOnly(False)
            self.register_ui_element(self.phase_plot_list, "plot_list")
        selected_labels_row.addWidget(self.phase_plot_list)
        select_labels_btn = QPushButton("Select Labels")
        select_labels_btn.clicked.connect(
            lambda: self.load_and_select_labels(self.phase_labels_file, self.phase_plot_list)
        )
        selected_labels_row.addWidget(select_labels_btn)
        plot_list_layout.addLayout(selected_labels_row)

        scroll_layout.addWidget(plot_list_group)

        # --- Color Dictionary ---
        color_dict_group = QGroupBox("Color Dictionary")
        color_dict_layout = QVBoxLayout()
        
        # Create or get existing color dict selector
        if not hasattr(self, 'plot_phase_col_dict'):
            self.plot_phase_col_dict = ColorDictSelector()
            self.register_ui_element(self.plot_phase_col_dict, "col_dict")
            
            # For HPA, pre-select Col_Dict_PMo if available
            if pom_type == "HPA":
                combo = self.plot_phase_col_dict.col_dict_combo
                for i in range(combo.count()):
                    if combo.itemText(i) == "Col_Dict_PMo":
                        combo.setCurrentIndex(i)
                        break
        
        color_dict_layout.addWidget(self.plot_phase_col_dict)
        color_dict_group.setLayout(color_dict_layout)
        scroll_layout.addWidget(color_dict_group)

        scroll_area.setWidget(scroll_widget)
        self.plot_phase_params_layout.addWidget(scroll_area)

        # Run / Stop buttons
        run_btn = QPushButton("▶ Generate Phase Plot")
        run_btn.setToolTip("Generate the Phase Diagram Plot")
        run_btn.setStyleSheet(ButtonStylesheets.RUN_BUTTON)
        run_btn.clicked.connect(self.run_plot_phase)

        stop_button = QPushButton("⏹ Stop")
        stop_button.setToolTip("Stop the current operation.")
        stop_button.setStyleSheet(ButtonStylesheets.STOP_BUTTON)
        stop_button.setEnabled(False)
        stop_button.clicked.connect(self.stop_function)

        self.plot_phase_run_btn = run_btn
        self.plot_phase_stop_btn = stop_button

        self.plot_phase_params_layout.addWidget(run_btn)
        self.plot_phase_params_layout.addWidget(stop_button)

    def run_plot_speciation(self):
        """Run the plot speciation diagram function"""
        self.clear_console()

        try:
            from utilities.plotting_gui import plot_speciation_run
            plot_function = plot_speciation_run
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import plotting function: {str(e)}")
            self.set_tab_error('plotting', 'plot_spec')
            return

        try:
            config_dict = self.get_gui_params()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to get GUI parameters: {str(e)}")
            self.set_tab_error('plotting', 'plot_spec')
            return

        # Validate required inputs
        if not config_dict['Preparation']['output_path'] or not config_dict['Preparation']['POM_system']:
            QMessageBox.critical(self, "Error", "Output path or POM system not selected.")
            self.set_tab_error('plotting', 'plot_spec')
            return

        try:
            self.update_console("Started plot speciation diagram", "plotting")
            self.set_tab_running('plotting', 'plot_spec')

            # Create and start function runner
            self.plot_spec_runner = POMSim_func_runner(plot_function, config_dict)
            self.plot_spec_runner.finished.connect(self.on_function_finished)
            self.plot_spec_runner.error.connect(self.on_function_error)
            self.plot_spec_runner.stopped.connect(self.on_function_stopped)
            # Connect console output to update console
            self.plot_spec_runner.console_output.connect(lambda text: self.update_console(text, "plotting"))
            
            self.plot_spec_run_btn.setEnabled(False)
            self.plot_spec_stop_btn.setEnabled(True)
            self.plot_spec_runner.start()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to start plotting function: {str(e)}")
            self.set_tab_error('plotting', 'plot_spec')
            self.plot_spec_run_btn.setEnabled(True)
            self.plot_spec_stop_btn.setEnabled(False)

    def run_plot_phase(self):
        """Run the plot phase diagram function"""
        self.clear_console()

        try:
            if self.pom_type == "IPA":
                from utilities.plotting_gui import plot_phase_ipa_run
                plot_function = plot_phase_ipa_run
            elif self.pom_type == "HPA":
                from utilities.plotting_gui import plot_phase_hpa_run
                plot_function = plot_phase_hpa_run
            else:
                QMessageBox.critical(self, "Error", "Invalid phase type selected.")
                self.set_tab_error('plotting', 'plot_phase')
                return
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to import plotting function: {str(e)}")
            self.set_tab_error('plotting', 'plot_phase')
            return

        try:
            config_dict = self.get_gui_params()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to get GUI parameters: {str(e)}")
            self.set_tab_error('plotting', 'plot_phase')
            return

        # Validate required inputs
        if not config_dict['Preparation']['output_path'] or not config_dict['Preparation']['POM_system']:
            QMessageBox.critical(self, "Error", "Output path or POM system not selected.")
            self.set_tab_error('plotting', 'plot_phase')
            return

        try:
            self.update_console("Started plot phase diagram", "plotting")
            self.set_tab_running('plotting', 'plot_phase')

            # Create and start function runner
            self.plot_phase_runner = POMSim_func_runner(plot_function, config_dict)
            self.plot_phase_runner.finished.connect(self.on_function_finished)
            self.plot_phase_runner.error.connect(self.on_function_error)
            self.plot_phase_runner.stopped.connect(self.on_function_stopped)
            # Connect console output to update console
            self.plot_phase_runner.console_output.connect(lambda text: self.update_console(text, "plotting"))
            
            self.plot_phase_run_btn.setEnabled(False)
            self.plot_phase_stop_btn.setEnabled(True)
            self.plot_phase_runner.start()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to start plotting function: {str(e)}")
            self.set_tab_error('plotting', 'plot_phase')
            self.plot_phase_run_btn.setEnabled(True)
            self.plot_phase_stop_btn.setEnabled(False)

class ButtonStylesheets:
    """
    Central repository for all button stylesheets used in the POM Simulator GUI.

    This class contains stylesheet constants for different button types and states,
    eliminating code duplication and making it easier to maintain consistent styling
    across the application.

    Implementation note: RUN_BUTTON and STOP_BUTTON share the same base appearance.
    The only difference is the disabled-state text color (black vs #666666).
    Both are derived from _BASE_ACTION_BUTTON to avoid repeating the common rules.
    """

    # Shared base for action buttons (Run and Stop).
    # Not intended for direct use; use RUN_BUTTON or STOP_BUTTON instead.
    _BASE_ACTION_BUTTON = """
        QPushButton {{
            font-size: 12pt;
            font-weight: bold;
            border: 2px solid #cccccc;
            background-color: #f4f4f4;
            color: black;
            border-radius: 5px;
            padding: 8px 16px;
            min-height: 20px;
        }}
        QPushButton:hover {{
            background-color: #f8f8f8;
            border-color: #333333;
        }}
        QPushButton:pressed {{
            background-color: #FFFFFF;
            border-color: #333333;
        }}
        QPushButton:disabled {{
            background-color: #cccccc;
            color: {disabled_color};
            border-color: #cccccc;
        }}
    """

    # Run button: disabled text is black (same as enabled) to keep it readable
    RUN_BUTTON = _BASE_ACTION_BUTTON.format(disabled_color="black")

    # Stop button: disabled text is dimmed to indicate the button is inactive
    STOP_BUTTON = _BASE_ACTION_BUTTON.format(disabled_color="#666666")

    # Browse button stylesheet (smaller, used for directory/file pickers)
    BROWSE_BUTTON = """
        QPushButton {
            font-size: 10pt;
            border: 1px solid #cccccc;
            background-color: #f4f4f4;
            color: black;
            border-radius: 3px;
            padding: 4px 8px;
        }
        QPushButton:hover {
            background-color: #f8f8f8;
            border-color: #333333;
        }
        QPushButton:pressed {
            background-color: #FFFFFF;
            border-color: #333333;
        }
        QPushButton:disabled {
            background-color: #cccccc;
            color: #666666;
            border-color: #cccccc;
        }
    """


class ColorDictSelector(QWidget):
    """
    Widget for selecting and visualizing color dictionaries from the database.

    Displays available color dictionaries (Col_Dict_PMo, Col_Dict_CU, etc.)
    and shows the label-to-color mappings with visual color swatches.
    """
    col_dict_changed = pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.available_col_dicts = {}
        self.current_col_dict = None
        self.current_col_dict_name = None
        self.init_ui()
        self.load_col_dicts_from_database()

    def init_ui(self):
        """Initialize the UI components"""
        layout = QVBoxLayout()

        # Selection dropdown
        selection_layout = QHBoxLayout()
        selection_layout.addWidget(QLabel("Color Dictionary:"))

        self.col_dict_combo = QComboBox()
        self.col_dict_combo.addItem("None (Default)")
        self.col_dict_combo.currentTextChanged.connect(self.on_col_dict_selected)
        selection_layout.addWidget(self.col_dict_combo)
        selection_layout.addStretch()

        layout.addLayout(selection_layout)

        # Color preview area
        preview_group = QGroupBox("Color Mappings")
        preview_layout = QVBoxLayout()

        # Scrollable area for color swatches
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setMinimumHeight(180)
        self.color_display_widget = QWidget()
        self.color_display_layout = QVBoxLayout(self.color_display_widget)
        scroll_area.setWidget(self.color_display_widget)
        preview_layout.addWidget(scroll_area)

        preview_group.setLayout(preview_layout)
        layout.addWidget(preview_group)

        self.setLayout(layout)

    def load_col_dicts_from_database(self):
        """Load all available color dictionaries from the database module"""
        try:
            from pomsimulator.modules.DataBase import color_dictionaries

            # Store available color dictionaries
            self.available_col_dicts = color_dictionaries

            # Add to combo box
            for col_dict_name in self.available_col_dicts.keys():
                self.col_dict_combo.addItem(col_dict_name)

        except ImportError as e:
            print(f"Warning: Could not import color dictionaries from database: {e}")
        except AttributeError as e:
            print(f"Warning: Color dictionary not found in database: {e}")

    def on_col_dict_selected(self, dict_name):
        """Handle color dictionary selection"""
        if dict_name == "None (Default)":
            self.current_col_dict = None
            self.current_col_dict_name = None
        else:
            self.current_col_dict = self.available_col_dicts.get(dict_name)
            self.current_col_dict_name = dict_name

        self.update_color_display()
        # Emit signal when color dictionary is updated
        self.col_dict_changed.emit(dict_name)

    def update_color_display(self):
        """Update the color swatch display in a grid layout"""
        # Clear existing display
        while self.color_display_layout.count():
            item = self.color_display_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

        if not self.current_col_dict:
            label = QLabel("No color dictionary selected")
            label.setStyleSheet("color: #999999; font-style: italic;")
            self.color_display_layout.addWidget(label)
            return

        # Create a grid layout for color swatches
        grid_layout = QGridLayout()
        grid_layout.setSpacing(15)
        grid_layout.setContentsMargins(10, 10, 10, 10)

        # Number of columns to display
        num_columns = 4
        row = 0
        col = 0

        # Display each color mapping in grid
        for label, color in self.current_col_dict.items():
            # Convert 8-character hex (with alpha) to 6-character hex for display
            # Format: #RRGGBBAA -> #RRGGBB
            if isinstance(color, str) and color.startswith('#'):
                if len(color) == 9:  # #RRGGBBAA format
                    display_color = color[:7]  # Take only #RRGGBB
                else:
                    display_color = color
            else:
                print(f"Warning: Invalid color format for label '{label}': {color}")
                continue

            # Create a container for each color item
            item_layout = QVBoxLayout()
            item_layout.setSpacing(5)
            item_layout.setContentsMargins(8, 8, 8, 8)

            # Color swatch
            swatch = QFrame()
            swatch.setFixedSize(100, 100)
            # Use the 6-character hex code for the stylesheet
            swatch.setStyleSheet(f"background-color: {display_color}; border: 2px solid #cccccc; border-radius: 5px;")
            item_layout.addWidget(swatch)

            # Label
            label_widget = QLabel(str(label))
            label_widget.setWordWrap(True)
            label_widget.setAlignment(Qt.AlignmentFlag.AlignCenter)
            label_widget.setStyleSheet("font-size: 10pt; font-weight: bold;")
            label_widget.setMaximumWidth(130)
            item_layout.addWidget(label_widget)

            # Color value (hex code) - show the original color code
            color_value = QLabel(str(color))
            color_value.setStyleSheet("color: #666666; font-family: monospace; font-size: 8pt;")
            color_value.setAlignment(Qt.AlignmentFlag.AlignCenter)
            color_value.setMaximumWidth(130)
            item_layout.addWidget(color_value)

            # Create container widget
            item_widget = QWidget()
            item_widget.setLayout(item_layout)
            item_widget.setStyleSheet(
                "border: 1px solid #e0e0e0; border-radius: 8px; padding: 5px; background-color: #f9f9f9;")
            item_widget.setMinimumWidth(150)
            # Add to grid
            grid_layout.addWidget(item_widget, row, col)

            # Move to next position
            col += 1
            if col >= num_columns:
                col = 0
                row += 1

        # Add the grid to the main layout
        grid_widget = QWidget()
        grid_widget.setLayout(grid_layout)

        # Add to scrollable area
        scroll_container = QWidget()
        scroll_container_layout = QVBoxLayout(scroll_container)
        scroll_container_layout.addWidget(grid_widget)
        scroll_container_layout.addStretch()

        self.color_display_layout.addWidget(scroll_container)

    def get_col_dict(self):
        """Get the current color dictionary for use in plotting functions"""
        return self.current_col_dict

    def get_col_dict_name(self):
        """Get the name of the current color dictionary"""
        return self.current_col_dict_name


class LabelSelectionDialog(QDialog):
    """Dialog for selecting labels from a labels file"""

    def __init__(self, labels_list, parent=None):
        super().__init__(parent)
        self.labels_list = labels_list
        self.selected_labels = []
        self.checkboxes = {}
        self.init_ui()

    def init_ui(self):
        """Initialize the dialog UI"""
        self.setWindowTitle("Select Labels")
        self.setGeometry(100, 100, 600, 500)

        layout = QVBoxLayout()

        # Checkboxes for labels
        checkbox_group = QGroupBox("Available Labels")
        checkbox_layout = QGridLayout()

        ncols = 4
        for idx, label in enumerate(sorted(self.labels_list)):
            row = idx // ncols
            col = idx % ncols

            checkbox = QCheckBox(label)
            checkbox.setChecked(False)
            self.checkboxes[label] = checkbox
            checkbox_layout.addWidget(checkbox, row, col)

        checkbox_group.setLayout(checkbox_layout)

        scroll_area = QScrollArea()
        scroll_area.setWidget(checkbox_group)
        scroll_area.setWidgetResizable(True)
        layout.addWidget(scroll_area)

        # Top buttons (Select All / Deselect All)
        top_buttons_layout = QHBoxLayout()

        select_all_btn = QPushButton("Select All")
        select_all_btn.clicked.connect(self.select_all)
        top_buttons_layout.addWidget(select_all_btn)

        deselect_all_btn = QPushButton("Deselect All")
        deselect_all_btn.clicked.connect(self.deselect_all)
        top_buttons_layout.addWidget(deselect_all_btn)

        layout.addLayout(top_buttons_layout)

        # Bottom buttons (Cancel / Include Selection)
        bottom_buttons_layout = QHBoxLayout()

        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        bottom_buttons_layout.addWidget(cancel_btn)

        include_btn = QPushButton("Include Selection")
        include_btn.clicked.connect(self.on_include_selection)
        bottom_buttons_layout.addWidget(include_btn)

        layout.addLayout(bottom_buttons_layout)

        self.setLayout(layout)

    def select_all(self):
        """Check all checkboxes"""
        for checkbox in self.checkboxes.values():
            checkbox.setChecked(True)

    def deselect_all(self):
        """Uncheck all checkboxes"""
        for checkbox in self.checkboxes.values():
            checkbox.setChecked(False)

    def get_selected_labels(self):
        """Get list of selected labels in sorted order"""
        selected = [label for label, checkbox in self.checkboxes.items() if checkbox.isChecked()]
        return sorted(selected)

    def on_include_selection(self):
        """Handle include selection button click with validation"""
        selected_labels = self.get_selected_labels()
        # Check if any labels are selected

        if not selected_labels:
            QMessageBox.warning(self, "Warning", "Please select at least one label.")
            return

        self.accept()


class CustomDoubleSpinBox(QDoubleSpinBox):
    """
    Custom QDoubleSpinBox that accepts both comma and dot as decimal separators.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Set locale that uses comma as decimal separator
        locale = QLocale()
        locale.setNumberOptions(QLocale.NumberOption.RejectGroupSeparator)
        self.setLocale(locale)

        # Enable keyboard tracking to respond to each keystroke
        self.setKeyboardTracking(True)

        # Allow intermediate editing states
        self.setStepType(QDoubleSpinBox.StepType.DefaultStepType)

    def validate(self, text, pos):
        # Handle empty text case
        if not text or text == '-':
            return (QValidator.Intermediate, text, pos)

        # Replace comma with dot for internal validation
        modified_text = text.replace(',', '.')

        # Try to convert to float to check if it's a valid number
        try:
            float(modified_text)
            return (QValidator.Acceptable, text, pos)
        except ValueError:
            # If it's not a valid float, check if it could be valid with more input
            if modified_text in ['-', '.', '-.', '-0.']:
                return (QValidator.Intermediate, text, pos)
            else:
                return (QValidator.Invalid, text, pos)

    def fixup(self, text):
        # If text is empty or invalid, return minimum value as string
        if not text:
            return str(self.minimum())

        # Replace comma with dot for fixup
        modified_text = text.replace(',', '.')

        # Try to convert to float
        try:
            value = float(modified_text)
            # Ensure value is within range
            value = max(min(value, self.maximum()), self.minimum())
            return str(value)
        except ValueError:
            return str(self.minimum())

    def valueFromText(self, text):
        # Handle empty text
        if not text:
            return self.minimum()

        # Replace comma with dot before converting to value
        modified_text = text.replace(',', '.')

        try:
            return float(modified_text)
        except ValueError:
            return self.minimum()

    def textFromValue(self, value):
        # Use the default text representation
        text = super().textFromValue(value)
        # If the locale uses comma, we keep it that way
        return text

    def stepBy(self, steps):
        # Ensure stepping works correctly
        super().stepBy(steps)

    def keyPressEvent(self, event):
        # Handle key press events
        super().keyPressEvent(event)


class ConsoleStream(StringIO):
    """Custom stream to capture print output and emit as signal"""

    def __init__(self, signal):
        super().__init__()
        self.signal = signal
        self._buffer = ""

    def write(self, text):
        if text:
            self._buffer += text
            # Only emit when we have a complete line (ends with newline)
            if '\n' in text:
                lines = self._buffer.split('\n')
                # Emit all complete lines
                for line in lines[:-1]:
                    if line.strip():  # Only emit non-empty lines
                        self.signal.emit(line)
                # Keep the last incomplete line in buffer
                self._buffer = lines[-1]
        return super().write(text)

    def flush(self):
        # Emit any remaining buffer content when flushed
        if self._buffer.strip():
            self.signal.emit(self._buffer)
            self._buffer = ""
        super().flush()


class POMSim_func_runner(QThread):
    """Thread class to run functions without blocking the GUI"""
    finished = pyqtSignal()
    error = pyqtSignal(str)
    console_output = pyqtSignal(str)  # Add signal for console output
    stopped = pyqtSignal()

    def __init__(self, function, config_dict):
        """
        Initialize the POMSim function runner thread.

        Creates a new thread instance for running POM Simulator functions asynchronously
        without blocking the GUI. Sets up a temporary stop file mechanism for graceful
        termination of long-running operations.

        Args:
            function (callable): The function to be executed in the separate thread.
                               This function should accept a configuration dictionary
                               as its parameter.
            config_dict (dict): Configuration dictionary containing parameters and
                              settings to be passed to the function during execution.

        Returns:
            None
        """
        super().__init__()
        self.function = function
        self.config_dict = config_dict

        self.is_running = False
        self.was_stopped = False
        self.stop_file = (tempfile.NamedTemporaryFile(delete=False))
        self.stop_file_name = self.stop_file.name
        self.stop_file.close()

    def run(self):
        """
        Execute the assigned function in a separate thread with output redirection.

        Runs the function specified during initialization while redirecting stdout and
        stderr to capture console output for GUI display. Implements a file-based stop
        mechanism and handles cleanup of resources. Emits appropriate signals for
        thread lifecycle events and error handling.

        Args:
            None

        Returns:
            None: This method doesn't return a value but emits signals to communicate
                  with the main thread (started, finished, error, console_output).
        """
        # Redirect stdout to capture prints
        original_stdout = sys.stdout
        original_stderr = sys.stderr

        console_stream = ConsoleStream(self.console_output)

        sys.stdout = console_stream
        sys.stderr = console_stream

        try:
            self.is_running = True
            self.was_stopped = False
            func = self.function

            # Clean up any existing stop file before starting
            try:
                if os.path.exists(self.stop_file_name):
                    os.unlink(self.stop_file_name)
                    # print(f"Cleaned up existing stop file: {self.stop_file_name}")
            except Exception as e:
                print(f"Warning: Could not clean up stop file: {e}")

            # Add stop file path to config dict so function can check if it should stop
            config_with_stop = self.config_dict.copy()
            config_with_stop['_stop_file'] = self.stop_file_name
            # print(f"Starting function: {self.function}")

            # Call the function with config_dict
            result = func(config_with_stop)

            # Check if we were stopped during execution
            if self.was_stopped or os.path.exists(self.stop_file_name):
                print("Function was stopped by user")
                self.stopped.emit()
                return

            if result == "Stopped":
                print(f"Function completed. Result: {result}")
                self.stopped.emit()
                return
            elif result == "Error" or result is None:
                print("Error running function")
                self.error.emit("Error running function")
                return
            else:
                print(f"Function completed. Result: {result}")
                self.finished.emit()

        except Exception as e:
            error_msg = f"Error running function: {str(e)}\n{traceback.format_exc()}"
            self.error.emit(error_msg)
            print(error_msg)
            return
        finally:
            # Restore original stdout
            sys.stdout = original_stdout
            sys.stderr = original_stderr

            self.is_running = False
            # Clean up the stop file
            try:
                if os.path.exists(self.stop_file_name):
                    os.unlink(self.stop_file_name)
            except Exception as e:
                print(f"Warning: Could not clean up stop file: {e}")

    def stop(self):
        """
        Signal the running function to stop execution using a file-based approach.

        Creates a stop file that the running function can check for to determine
        if it should terminate gracefully. This provides a mechanism for stopping
        long-running operations without forcefully terminating the thread.

        Args:
            None

        Returns:
            None: This method doesn't return a value but creates a stop file
                  and updates the is_running flag to signal termination.
        """
        self.is_running = False
        self.was_stopped = True
        # Create the stop file to signal stopping
        try:
            with open(self.stop_file_name, 'w') as f:
                f.write('stop')
            print(f"Stop signal written to: {self.stop_file_name}")
        except Exception as e:
            print(f"Error writing stop file: {e}")


def main():
    """
    Initialize and run the POM Simulator GUI application.

    This function serves as the entry point for the POM Simulator GUI application.
    It creates a QApplication instance, initializes the main window of the
    application (POMSimulatorGUI), displays it, and starts the application's
    event loop. The function will only return when the application is closed,
    at which point it ensures proper termination with the appropriate exit code.

    Parameters:
        None

    Returns:
        None: This function doesn't return as it calls sys.exit() to terminate
              the program with the exit code from app.exec_().
    """
    app = QApplication(sys.argv)

    # Set application style for better appearance
    app.setStyle('Fusion')

    # Create and show the main window
    window = POMSimulatorGUI()
    window.show()

    # Process events to ensure the window is displayed
    app.processEvents()

    # Start the event loop
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()











    sys.exit(app.exec_())


if __name__ == "__main__":
    main()












