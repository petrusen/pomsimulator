"""
DPI Scaling Utilities for POMSimulator

This module provides utilities for making the PyQt5 GUI DPI-aware and visually
consistent across different platforms and display settings.

The scaling system uses the monitor's logical DPI to compute a scale factor
relative to a reference DPI of 96 (standard Windows DPI). This ensures that
the GUI maintains approximately the same physical appearance regardless of
screen resolution, size, DPI settings, or operating system.
"""

from PyQt5.QtWidgets import QApplication
from PyQt5.QtGui import QFont
from PyQt5.QtCore import QSize, Qt


class DPIScaler:
    """
    Centralized DPI scaling utility class.
    
    This class computes and applies DPI-aware scaling throughout the application.
    It uses Qt's logical DPI information to maintain consistent physical sizing
    across different displays and platforms.
    """
    
    _instance = None
    _scale_factor = None
    
    def __new__(cls):
        """Singleton pattern to ensure consistent scaling throughout the app."""
        if cls._instance is None:
            cls._instance = super(DPIScaler, cls).__new__(cls)
        return cls._instance
    
    def __init__(self):
        """Initialize the DPI scaler if not already initialized."""
        if self._scale_factor is None:
            self._compute_scale_factor()
    
    def _compute_scale_factor(self):
        """
        Compute the DPI scale factor based on the primary screen's logical DPI.
        
        Uses 96 DPI as the reference value (standard Windows DPI).
        The scale factor is clamped to a reasonable range (0.9-2.0) to avoid
        extreme cases that could make the UI unusable.
        """
        app = QApplication.instance()
        if app is None:
            # Fallback if no QApplication exists yet
            self._scale_factor = 1.0
            return
        
        try:
            # Get the primary screen
            primary_screen = app.primaryScreen()
            if primary_screen is None:
                self._scale_factor = 1.0
                return
            
            # Get logical DPI (accounts for system scaling settings)
            logical_dpi = primary_screen.logicalDotsPerInch()
            
            # Compute scale factor relative to 96 DPI reference
            raw_scale = logical_dpi / 96.0
            
            # Clamp to reasonable range to avoid extreme cases
            self._scale_factor = max(0.9, min(2.0, raw_scale))
            
        except Exception:
            # Fallback to no scaling if anything goes wrong
            self._scale_factor = 1.0
    
    def get_scale_factor(self):
        """
        Get the current DPI scale factor.
        
        Returns:
            float: The scale factor to apply to pixel values.
        """
        if self._scale_factor is None:
            self._compute_scale_factor()
        return self._scale_factor
    
    def scale(self, value):
        """
        Scale a single pixel value.
        
        Args:
            value (int or float): The pixel value to scale.
            
        Returns:
            int: The scaled pixel value, rounded to nearest integer.
        """
        if value == 0:
            return 0
        return int(round(value * self.get_scale_factor()))
    
    def scale_size(self, width, height):
        """
        Scale width and height values.
        
        Args:
            width (int): The width in pixels to scale.
            height (int): The height in pixels to scale.
            
        Returns:
            tuple: (scaled_width, scaled_height) as integers.
        """
        return (self.scale(width), self.scale(height))
    
    def scale_qsize(self, width, height):
        """
        Create a scaled QSize object.
        
        Args:
            width (int): The width in pixels to scale.
            height (int): The height in pixels to scale.
            
        Returns:
            QSize: A QSize object with scaled dimensions.
        """
        return QSize(self.scale(width), self.scale(height))
    
    def scale_font(self, base_points):
        """
        Create a scaled font with point-based sizing.
        
        Args:
            base_points (int or float): The base font size in points.
            
        Returns:
            QFont: A QFont object with scaled point size.
        """
        font = QFont()
        # Scale font size but ensure it's at least 6 points
        scaled_points = max(6, int(round(base_points * self.get_scale_factor())))
        font.setPointSize(scaled_points)
        return font
    
    def scale_monospace_font(self, base_points):
        """
        Create a scaled monospace font with point-based sizing.
        
        Args:
            base_points (int or float): The base font size in points.
            
        Returns:
            QFont: A monospace QFont object with scaled point size.
        """
        font = QFont("Monospace")
        font.setStyleHint(QFont.Monospace)
        # Scale font size but ensure it's at least 6 points
        scaled_points = max(6, int(round(base_points * self.get_scale_factor())))
        font.setPointSize(scaled_points)
        return font


# Global instance for easy access
_dpi_scaler = DPIScaler()


def scale(value):
    """
    Scale a single pixel value using the global DPI scaler.
    
    Args:
        value (int or float): The pixel value to scale.
        
    Returns:
        int: The scaled pixel value.
    """
    return _dpi_scaler.scale(value)


def scale_size(width, height):
    """
    Scale width and height values using the global DPI scaler.
    
    Args:
        width (int): The width in pixels to scale.
        height (int): The height in pixels to scale.
        
    Returns:
        tuple: (scaled_width, scaled_height) as integers.
    """
    return _dpi_scaler.scale_size(width, height)


def scale_qsize(width, height):
    """
    Create a scaled QSize object using the global DPI scaler.
    
    Args:
        width (int): The width in pixels to scale.
        height (int): The height in pixels to scale.
        
    Returns:
        QSize: A QSize object with scaled dimensions.
    """
    return _dpi_scaler.scale_qsize(width, height)


def scale_font(base_points):
    """
    Create a scaled font using the global DPI scaler.
    
    Args:
        base_points (int or float): The base font size in points.
        
    Returns:
        QFont: A QFont object with scaled point size.
    """
    return _dpi_scaler.scale_font(base_points)


def scale_monospace_font(base_points):
    """
    Create a scaled monospace font using the global DPI scaler.
    
    Args:
        base_points (int or float): The base font size in points.
        
    Returns:
        QFont: A monospace QFont object with scaled point size.
    """
    return _dpi_scaler.scale_monospace_font(base_points)


def get_scale_factor():
    """
    Get the current DPI scale factor.
    
    Returns:
        float: The scale factor being applied to pixel values.
    """
    return _dpi_scaler.get_scale_factor()


def scale_css_font_size(css_text):
    """
    Scale font sizes in CSS text.
    
    This function finds font-size declarations in CSS and scales them
    according to the current DPI scale factor.
    
    Args:
        css_text (str): CSS text containing font-size declarations.
        
    Returns:
        str: CSS text with scaled font sizes.
    """
    import re
    
    def replace_font_size(match):
        size_str = match.group(1)
        try:
            size = float(size_str)
            scaled_size = max(6, int(round(size * get_scale_factor())))
            return f"font-size: {scaled_size}pt"
        except ValueError:
            return match.group(0)  # Return original if parsing fails
    
    # Pattern to match font-size: XXpt declarations
    pattern = r'font-size:\s*(\d+(?:\.\d+)?)pt'
    return re.sub(pattern, replace_font_size, css_text)


def enable_high_dpi_support():
    """
    Enable Qt High-DPI support.
    
    This function should be called BEFORE creating the QApplication instance.
    It enables Qt's built-in high-DPI scaling mechanisms for better cross-platform
    compatibility.
    
    Note: This function has no effect if called after QApplication is created.
    """
    try:
        # Enable High DPI display with PyQt5
        QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
        QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
        
        # For PyQt5 5.14+, also set the high DPI scale factor rounding policy
        try:
            if hasattr(Qt, 'HighDpiScaleFactorRoundingPolicy'):
                QApplication.setHighDpiScaleFactorRoundingPolicy(
                    Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
                )
        except AttributeError:
            # Older PyQt5 versions don't have this attribute
            pass
            
    except Exception:
        # If anything fails, continue without high DPI support
        # This ensures compatibility with older PyQt5 versions
        pass