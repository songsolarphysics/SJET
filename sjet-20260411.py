import numpy as np
import matplotlib

import matplotlib.pyplot as plt
from scipy import ndimage
from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from skimage import filters, morphology, feature, measure
from skimage.draw import polygon as draw_polygon
import cv2
import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import RectangleSelector, PolygonSelector
import sunpy.map
import astropy.units as u
from matplotlib.colors import LogNorm

class InteractiveThresholdingApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Solar Jet Extraction Tool")
        self.root.geometry("1200x800")

        # Data storage
        self.fits_file = None
        self.data = None
        self.header = None
        self.sunpy_map = None
        self.normalized_data = None
        self.binary_mask = None
        self.Extracted_jet = None
        self.edges = None
        self._data_cache = None  # Cache for processing
        self.max_region_only = tk.BooleanVar(value=False)  # Variable for max region option

        # Region merging variables
        self.merge_regions = tk.BooleanVar(value=False)  # Flag for region merging
        self.merge_distance = tk.IntVar(value=10)  # Distance threshold for merging
        self.merge_size_threshold = tk.IntVar(value=50)  # Minimum size for regions to be considered
        self.force_merge_all = tk.BooleanVar(value=False)  # New flag for force merging all regions
        self.merge_strength = tk.DoubleVar(value=0.15)  # New variable for merge connection strength

        # ROI selection variables
        self.roi_active = tk.BooleanVar(value=False)  # Flag for ROI selection mode
        self.roi_mode = tk.StringVar(value="rectangle")  # ROI selection mode: rectangle or polygon
        self.roi_coords = None  # Will store the coordinates of selected ROI
        self.roi_polygon_points = []  # Store polygon points for polygon mode
        self.roi_mask = None    # Will store the mask for ROI
        self.rect_selector = None  # Will store the RectangleSelector object
        self.polygon_selector = None  # Will store the PolygonSelector object
        self.roi_fig = None  # Store the ROI selection figure
        self.roi_ax = None  # Store the ROI selection axes

        # Visualization caching
        self.norm_log = None
        self.cbar1 = None
        self.cbar3 = None

        # Create UI
        self.create_ui()

    def create_ui(self):
        # Main frame
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Create a scrollable frame for controls
        control_container = ttk.Frame(main_frame)
        control_container.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)

        # Create canvas and scrollbar
        control_canvas = tk.Canvas(control_container)
        control_scrollbar = ttk.Scrollbar(control_container, orient="vertical", command=control_canvas.yview)

        # Configure the canvas
        control_canvas.configure(yscrollcommand=control_scrollbar.set)

        # Pack canvas and scrollbar
        control_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        control_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Create the scrollable frame for controls
        control_frame = ttk.LabelFrame(control_canvas, text="Controls")

        # Create window in canvas
        control_frame_window = control_canvas.create_window((0, 0), window=control_frame, anchor="nw")

        # Configure the interior frame to update the scrollregion
        def _configure_interior(event):
            size = (control_frame.winfo_reqwidth(), control_frame.winfo_reqheight())
            control_canvas.config(scrollregion=control_canvas.bbox("all"))
            if control_frame.winfo_reqwidth() != control_canvas.winfo_width():
                control_canvas.config(width=control_frame.winfo_reqwidth())

        def _configure_canvas(event):
            if control_frame.winfo_reqwidth() != control_canvas.winfo_width():
                control_canvas.itemconfigure(control_frame_window, width=control_canvas.winfo_width())

        # Bind events
        control_frame.bind('<Configure>', _configure_interior)
        control_canvas.bind('<Configure>', _configure_canvas)    
        
        # File selection
        file_frame = ttk.Frame(control_frame)
        file_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(file_frame, text="Open FITS File", command=self.open_fits_file).pack(fill=tk.X)
        ttk.Label(file_frame, text="Current File:").pack(anchor=tk.W, pady=(10, 0))
        self.file_label = ttk.Label(file_frame, text="None selected")
        self.file_label.pack(anchor=tk.W)

        # ROI Selection Frame
        roi_frame = ttk.LabelFrame(control_frame, text="Region of Interest (ROI)")
        roi_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Checkbutton(roi_frame, text="Enable ROI Selection",
                       variable=self.roi_active,
                       command=self.toggle_roi_selection).pack(anchor=tk.W, padx=5, pady=5)

        # ROI Mode Selection
        roi_mode_frame = ttk.Frame(roi_frame)
        roi_mode_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(roi_mode_frame, text="ROI Mode:").pack(side=tk.LEFT, padx=5)
        ttk.Radiobutton(roi_mode_frame, text="Rectangle", variable=self.roi_mode,
                        value="rectangle").pack(side=tk.LEFT, padx=2)
        ttk.Radiobutton(roi_mode_frame, text="Polygon", variable=self.roi_mode,
                        value="polygon").pack(side=tk.LEFT, padx=2)

        roi_buttons_frame = ttk.Frame(roi_frame)
        roi_buttons_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Button(roi_buttons_frame, text="Select ROI",
                  command=self.start_roi_selection).pack(side=tk.LEFT, padx=2)
        ttk.Button(roi_buttons_frame, text="Clear ROI",
                  command=self.clear_roi).pack(side=tk.LEFT, padx=2)
        ttk.Button(roi_buttons_frame, text="Apply ROI",
                  command=self.apply_roi).pack(side=tk.LEFT, padx=2)

        self.roi_status_label = ttk.Label(roi_frame, text="ROI Status: No ROI selected")
        self.roi_status_label.pack(anchor=tk.W, padx=5, pady=5)

        # Thresholding method
        method_frame = ttk.LabelFrame(control_frame, text="Thresholding Method")
        method_frame.pack(fill=tk.X, padx=5, pady=5)

        self.method_var = tk.StringVar(value="manual")
        ttk.Radiobutton(method_frame, text="Manual Threshold", variable=self.method_var,
                        value="manual", command=self.update_method).pack(anchor=tk.W)
        ttk.Radiobutton(method_frame, text="Otsu's Method", variable=self.method_var,
                        value="otsu", command=self.update_method).pack(anchor=tk.W)
        ttk.Radiobutton(method_frame, text="Adaptive Threshold", variable=self.method_var,
                        value="adaptive", command=self.update_method).pack(anchor=tk.W)
        ttk.Radiobutton(method_frame, text="Percentile", variable=self.method_var,
                        value="percentile", command=self.update_method).pack(anchor=tk.W)
        ttk.Radiobutton(method_frame, text="Log Extracted", variable=self.method_var,
                        value="log", command=self.update_method).pack(anchor=tk.W)

        # Manual threshold slider
        threshold_frame = ttk.LabelFrame(control_frame, text="Manual Threshold")
        threshold_frame.pack(fill=tk.X, padx=5, pady=5)

        self.threshold_var = tk.DoubleVar(value=50)
        self.threshold_slider = ttk.Scale(threshold_frame, from_=0, to=100,
                                         variable=self.threshold_var, orient=tk.HORIZONTAL,
                                         command=self.update_threshold)
        self.threshold_slider.pack(fill=tk.X, padx=5, pady=5)

        self.threshold_label = ttk.Label(threshold_frame, text="Threshold: 50%")
        self.threshold_label.pack(anchor=tk.W, padx=5)

        # Adaptive threshold parameters
        adaptive_frame = ttk.LabelFrame(control_frame, text="Adaptive Parameters")
        adaptive_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Label(adaptive_frame, text="Block Size:").pack(anchor=tk.W, padx=5)
        self.block_size_var = tk.IntVar(value=51)
        self.block_size_slider = ttk.Scale(adaptive_frame, from_=3, to=101,
                                          variable=self.block_size_var, orient=tk.HORIZONTAL,
                                          command=self.update_block_size)
        self.block_size_slider.pack(fill=tk.X, padx=5, pady=5)
        self.block_size_label = ttk.Label(adaptive_frame, text="Block Size: 51")
        self.block_size_label.pack(anchor=tk.W, padx=5)

        ttk.Label(adaptive_frame, text="C Value:").pack(anchor=tk.W, padx=5)
        self.c_value_var = tk.IntVar(value=2)
        self.c_value_slider = ttk.Scale(adaptive_frame, from_=0, to=20,
                                       variable=self.c_value_var, orient=tk.HORIZONTAL,
                                       command=self.update_c_value)
        self.c_value_slider.pack(fill=tk.X, padx=5, pady=5)
        self.c_value_label = ttk.Label(adaptive_frame, text="C Value: 2")
        self.c_value_label.pack(anchor=tk.W, padx=5)

        # Percentile threshold
        percentile_frame = ttk.LabelFrame(control_frame, text="Percentile")
        percentile_frame.pack(fill=tk.X, padx=5, pady=5)

        self.percentile_var = tk.DoubleVar(value=95)
        self.percentile_slider = ttk.Scale(percentile_frame, from_=50, to=99.9,
                                          variable=self.percentile_var, orient=tk.HORIZONTAL,
                                          command=self.update_percentile)
        self.percentile_slider.pack(fill=tk.X, padx=5, pady=5)
        self.percentile_label = ttk.Label(percentile_frame, text="Percentile: 95%")
        self.percentile_label.pack(anchor=tk.W, padx=5)

        # Log Extracted parameters
        log_frame = ttk.LabelFrame(control_frame, text="Log Extracted")
        log_frame.pack(fill=tk.X, padx=5, pady=5)

        self.log_factor_var = tk.DoubleVar(value=1.5)
        self.log_factor_slider = ttk.Scale(log_frame, from_=0.1, to=5.0,
                                          variable=self.log_factor_var, orient=tk.HORIZONTAL,
                                          command=self.update_log_factor)
        self.log_factor_slider.pack(fill=tk.X, padx=5, pady=5)
        self.log_factor_label = ttk.Label(log_frame, text="Factor: 1.5")
        self.log_factor_label.pack(anchor=tk.W, padx=5)

        # Morphological operations
        morph_frame = ttk.LabelFrame(control_frame, text="Morphological Operations")
        morph_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Label(morph_frame, text="Opening Size:").pack(anchor=tk.W, padx=5)
        self.opening_var = tk.IntVar(value=2)
        self.opening_slider = ttk.Scale(morph_frame, from_=0, to=10,
                                       variable=self.opening_var, orient=tk.HORIZONTAL,
                                       command=self.update_opening)
        self.opening_slider.pack(fill=tk.X, padx=5, pady=5)

        ttk.Label(morph_frame, text="Closing Size:").pack(anchor=tk.W, padx=5)
        self.closing_var = tk.IntVar(value=3)
        self.closing_slider = ttk.Scale(morph_frame, from_=0, to=10,
                                       variable=self.closing_var, orient=tk.HORIZONTAL,
                                       command=self.update_closing)
        self.closing_slider.pack(fill=tk.X, padx=5, pady=5)

        # Maximum region option
        max_region_frame = ttk.LabelFrame(control_frame, text="Maximum Region Filter")
        max_region_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Checkbutton(max_region_frame, text="Keep only region with maximum intensity",
                       variable=self.max_region_only,
                       command=self.update_max_region).pack(anchor=tk.W, padx=5, pady=5)

        # Intensity percentile for max region
        ttk.Label(max_region_frame, text="Max Region Intensity Threshold:").pack(anchor=tk.W, padx=5)
        self.max_intensity_var = tk.DoubleVar(value=90)
        self.max_intensity_slider = ttk.Scale(max_region_frame, from_=50, to=99.9,
                                            variable=self.max_intensity_var, orient=tk.HORIZONTAL,
                                            command=self.update_max_intensity)
        self.max_intensity_slider.pack(fill=tk.X, padx=5, pady=5)
        self.max_intensity_label = ttk.Label(max_region_frame, text="Intensity Threshold: 90%")
        self.max_intensity_label.pack(anchor=tk.W, padx=5)

        # Region Merging Frame
        merge_frame = ttk.LabelFrame(control_frame, text="Region Merging")
        merge_frame.pack(fill=tk.X, padx=5, pady=5)

        ttk.Checkbutton(merge_frame, text="Enable automatic region merging",
                       variable=self.merge_regions,
                       command=self.update_merge_regions).pack(anchor=tk.W, padx=5, pady=5)

        # Force merge all regions option
        ttk.Checkbutton(merge_frame, text="Force merge all regions (ignore distance)",
                       variable=self.force_merge_all,
                       command=self.update_force_merge).pack(anchor=tk.W, padx=5, pady=5)

        # Merge connection strength slider - NEW
        ttk.Label(merge_frame, text="Merge Connection Strength:").pack(anchor=tk.W, padx=5)
        self.merge_strength_slider = ttk.Scale(merge_frame, from_=0.01, to=0.3,
                                             variable=self.merge_strength, orient=tk.HORIZONTAL,
                                             command=self.update_merge_strength)
        self.merge_strength_slider.pack(fill=tk.X, padx=5, pady=5)
        self.merge_strength_label = ttk.Label(merge_frame, text="Connection Strength: 0.15 (Use Convex Hull)")
        self.merge_strength_label.pack(anchor=tk.W, padx=5)

        # Merge distance slider (only used when force merge is off)
        ttk.Label(merge_frame, text="Merge Distance Threshold:").pack(anchor=tk.W, padx=5)
        self.merge_distance_slider = ttk.Scale(merge_frame, from_=1, to=50,
                                            variable=self.merge_distance, orient=tk.HORIZONTAL,
                                            command=self.update_merge_distance)
        self.merge_distance_slider.pack(fill=tk.X, padx=5, pady=5)
        self.merge_distance_label = ttk.Label(merge_frame, text="Distance Threshold: 10 pixels")
        self.merge_distance_label.pack(anchor=tk.W, padx=5)

        # Merge size threshold slider
        ttk.Label(merge_frame, text="Minimum Region Size:").pack(anchor=tk.W, padx=5)
        self.merge_size_slider = ttk.Scale(merge_frame, from_=1, to=500,
                                         variable=self.merge_size_threshold, orient=tk.HORIZONTAL,
                                         command=self.update_merge_size)
        self.merge_size_slider.pack(fill=tk.X, padx=5, pady=5)
        self.merge_size_label = ttk.Label(merge_frame, text="Minimum Size: 50 pixels")
        self.merge_size_label.pack(anchor=tk.W, padx=5)

        # Save button
        save_frame = ttk.Frame(control_frame)
        save_frame.pack(fill=tk.X, padx=5, pady=15)
        ttk.Button(save_frame, text="Save Results", command=self.save_results).pack(fill=tk.X)

        # Right panel for visualization
        viz_frame = ttk.LabelFrame(main_frame, text="Visualization")
        viz_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Create matplotlib figure
        self.fig, self.axes = plt.subplots(2, 2, figsize=(8, 6))
        self.fig.tight_layout(pad=3.0)

        # Embed matplotlib figure in tkinter
        self.canvas = FigureCanvasTkAgg(self.fig, master=viz_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Initialize with empty plots
        for ax in self.axes.flat:
            ax.set_xticks([])
            ax.set_yticks([])

        self.axes[0, 0].set_title("Original Data")
        self.axes[0, 1].set_title("Binary Mask")
        self.axes[1, 0].set_title("Extracted Jet")
        self.axes[1, 1].set_title("Edge Detection")

        # Update UI based on initial method
        self.update_method()
        self.update_force_merge()

    # Merge strength control method
    def update_merge_strength(self, *args):
        value = self.merge_strength.get()
        if value > 0.1:
            self.merge_strength_label.config(text=f"Connection Strength: {value:.2f} (Use Convex Hull)")
        else:
            self.merge_strength_label.config(text=f"Connection Strength: {value:.2f} (No Convex Hull)")
            
        if self.data is not None and self.merge_regions.get() and self.force_merge_all.get():
            self.process_image()

    # Force merge update method
    def update_force_merge(self, *args):
        if self.force_merge_all.get():
            self.merge_strength_slider.state(["!disabled"])
            self.merge_strength_label.state(["!disabled"])
            self.merge_distance_slider.state(["disabled"])
            self.merge_distance_label.config(text="Distance Threshold: Not Used")
        else:
            self.merge_strength_slider.state(["disabled"])
            self.merge_strength_label.state(["disabled"])
            self.merge_distance_slider.state(["!disabled"])
            self.merge_distance_label.config(text=f"Distance Threshold: {self.merge_distance.get()} pixels")

        if self.data is not None and self.merge_regions.get():
            self.process_image()

    # Region merging control methods
    def update_merge_regions(self, *args):
        if self.data is not None:
            self.process_image()

    def update_merge_distance(self, *args):
        value = self.merge_distance.get()
        if not self.force_merge_all.get():
            self.merge_distance_label.config(text=f"Distance Threshold: {value} pixels")
        if self.data is not None and self.merge_regions.get() and not self.force_merge_all.get():
            self.process_image()

    def update_merge_size(self, *args):
        value = self.merge_size_threshold.get()
        self.merge_size_label.config(text=f"Minimum Size: {value} pixels")
        if self.data is not None and self.merge_regions.get():
            self.process_image()

    def open_fits_file(self):
        file_path = filedialog.askopenfilename(
            title="Select FITS File",
            filetypes=[("FITS files", "*.fits"), ("All files", "*.*")]
        )

        if file_path:
            try:
                # Use SunPy to read FITS file
                self.sunpy_map = sunpy.map.Map(file_path)
                self.data = self.sunpy_map.data
                self.header = self.sunpy_map.meta

                # Handle NaN values
                self.data = np.nan_to_num(self.data)

                # Normalize data for display
                self.normalized_data = ((self.data - np.min(self.data)) /
                                       (np.max(self.data) - np.min(self.data)) * 255).astype(np.uint8)

                # Initialize ROI mask as full image
                self.roi_mask = np.ones_like(self.data, dtype=bool)

                # Initialize data cache
                self._data_cache = self.data.copy()

                # Update file label
                self.file_label.config(text=os.path.basename(file_path))
                self.fits_file = file_path

                # Process image with current settings
                self.process_image()

            except Exception as e:
                messagebox.showerror("Error", f"Failed to open FITS file: {str(e)}")

    # ROI SELECTION METHODS
    def toggle_roi_selection(self):
        """Toggle ROI selection mode"""
        if self.data is None:
            messagebox.showwarning("Warning", "Please load a FITS file first")
            self.roi_active.set(False)
            return

        if self.roi_active.get():
            self.roi_status_label.config(text="ROI Status: Selection mode active")
        else:
            self.roi_status_label.config(text="ROI Status: Selection mode inactive")

    def start_roi_selection(self):
        if self.data is None:
            messagebox.showwarning("Warning", "Please load a FITS file first")
            return
        self.roi_active.set(True)
        # Close previous ROI window
        if hasattr(self, 'roi_win') and self.roi_win is not None:
            try:
                self.roi_win.destroy()
            except:
                pass
        # Create a Toplevel window sharing the same tkinter mainloop
        self.roi_win = tk.Toplevel(self.root)
        self.roi_win.title("ROI Selection")
        self.roi_win.geometry("900x750")
        # Embed matplotlib figure in tkinter
        self.roi_fig, self.roi_ax = plt.subplots(figsize=(9, 7))
        self.roi_ax.imshow(self.data, cmap='inferno', origin='lower',
                           norm=LogNorm(vmin=max(1, np.nanmin(self.data[self.data > 0])),
                                       vmax=np.nanmax(self.data)))
        roi_canvas = FigureCanvasTkAgg(self.roi_fig, master=self.roi_win)
        roi_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        # Add navigation toolbar for zoom/pan
        from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
        toolbar = NavigationToolbar2Tk(roi_canvas, self.roi_win)
        toolbar.update()
        mode = self.roi_mode.get()
        if mode == "rectangle":
            self.roi_ax.set_title(
                "Drag to select a rectangular ROI, then click Confirm", fontsize=11)
            self.rect_selector = RectangleSelector(
                self.roi_ax, self.roi_callback,
                useblit=False,
                button=[1],
                minspanx=5, minspany=5,
                spancoords='pixels',
                interactive=True,
                props=dict(facecolor='cyan', edgecolor='lime',
                           alpha=0.2, fill=True)
            )
        else:
            self.roi_ax.set_title(
                "Click to add polygon vertices, press Enter to confirm, Esc to cancel", fontsize=11)
            self.roi_polygon_points = []
            self.polygon_selector = PolygonSelector(
                self.roi_ax, self.polygon_callback,
                useblit=False,
                props=dict(color='lime', linewidth=2)
            )
        roi_canvas.draw()
        # Confirm button: close window and update status
        def on_confirm():
            if mode == "rectangle" and self.roi_coords is not None:
                self.roi_status_label.config(
                    text=f"ROI Status: Rectangle selected {self.roi_coords} — click Apply ROI to apply")
            elif mode == "polygon" and len(self.roi_polygon_points) >= 3:
                self.roi_status_label.config(
                    text=f"ROI Status: Polygon selected with {len(self.roi_polygon_points)} vertices — click Apply ROI to apply")
            else:
                self.roi_status_label.config(text="ROI Status: No valid selection detected")
            self.roi_win.destroy()
            self.roi_win = None
            plt.close(self.roi_fig)
        tk.Button(
            self.roi_win, text="✓  Confirm ROI",
            bg='#226622', fg='white', font=('Arial', 12, 'bold'),
            command=on_confirm
        ).pack(pady=8)    
    
    
    
    def roi_callback(self, eclick, erelease):
        """Callback function for rectangle ROI selection"""
        x1, y1 = int(eclick.xdata), int(eclick.ydata)
        x2, y2 = int(erelease.xdata), int(erelease.ydata)

        # Ensure coordinates are within image bounds
        height, width = self.data.shape
        x1 = max(0, min(x1, width-1))
        y1 = max(0, min(y1, height-1))
        x2 = max(0, min(x2, width-1))
        y2 = max(0, min(y2, height-1))

        # Store the ROI coordinates
        self.roi_coords = (min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2))

        # Update status label
        self.roi_status_label.config(
            text=f"ROI Status: Selected ({self.roi_coords[0]},{self.roi_coords[1]}) to ({self.roi_coords[2]},{self.roi_coords[3]})"
        )

    def polygon_callback(self, vertices):
        """Callback function for polygon ROI selection"""
        if len(vertices) < 3:
            messagebox.showwarning("Warning", "Polygon must have at least 3 vertices")
            return
        
        # Store polygon vertices
        self.roi_polygon_points = vertices
        
        # Update status label
        self.roi_status_label.config(
            text=f"ROI Status: Polygon selected with {len(vertices)} vertices"
        )

    def clear_roi(self):
        """Clear the current ROI selection"""
        if self.data is None:
            return

        self.roi_coords = None
        self.roi_polygon_points = []
        self.roi_mask = np.ones_like(self.data, dtype=bool)
        self.roi_status_label.config(text="ROI Status: No ROI selected")

        # Reprocess the image without ROI constraint
        self.process_image()

    def apply_roi(self):
        """Apply the selected ROI to the processing"""
        mode = self.roi_mode.get()
        
        if mode == "rectangle":
            if self.roi_coords is None:
                messagebox.showwarning("Warning", "Please select a ROI first")
                return

            # Create ROI mask efficiently
            self.roi_mask = np.zeros_like(self.data, dtype=bool)
            x1, y1, x2, y2 = self.roi_coords
            self.roi_mask[y1:y2+1, x1:x2+1] = True

            # Update status
            self.roi_status_label.config(text=f"ROI Status: Applied ({x1},{y1}) to ({x2},{y2})")
            
        else:  # polygon mode
            if len(self.roi_polygon_points) < 3:
                messagebox.showwarning("Warning", "Please select a polygon ROI with at least 3 points")
                return

            # Create polygon mask
            self.roi_mask = np.zeros_like(self.data, dtype=bool)
            vertices = np.array(self.roi_polygon_points)
            
            # Use skimage's polygon function to create mask
            rr, cc = draw_polygon(vertices[:, 1], vertices[:, 0], self.data.shape)
            self.roi_mask[rr, cc] = True

            # Update status
            self.roi_status_label.config(text=f"ROI Status: Polygon applied with {len(self.roi_polygon_points)} vertices")

        # Process the image
        self.process_image(roi_changed=True)

    def update_method(self):
        method = self.method_var.get()

        # Enable/disable relevant controls based on selected method
        if method == "manual":
            self.threshold_slider.state(["!disabled"])
            self.threshold_label.state(["!disabled"])
        else:
            self.threshold_slider.state(["disabled"])
            self.threshold_label.state(["disabled"])

        if method == "adaptive":
            self.block_size_slider.state(["!disabled"])
            self.block_size_label.state(["!disabled"])
            self.c_value_slider.state(["!disabled"])
            self.c_value_label.state(["!disabled"])
        else:
            self.block_size_slider.state(["disabled"])
            self.block_size_label.state(["disabled"])
            self.c_value_slider.state(["disabled"])
            self.c_value_label.state(["disabled"])

        if method == "percentile":
            self.percentile_slider.state(["!disabled"])
            self.percentile_label.state(["!disabled"])
        else:
            self.percentile_slider.state(["disabled"])
            self.percentile_label.state(["disabled"])

        if method == "log":
            self.log_factor_slider.state(["!disabled"])
            self.log_factor_label.state(["!disabled"])
        else:
            self.log_factor_slider.state(["disabled"])
            self.log_factor_label.state(["disabled"])

        # Process image with new method
        if self.data is not None:
            self.process_image()

    def update_threshold(self, *args):
        value = self.threshold_var.get()
        self.threshold_label.config(text=f"Threshold: {value:.1f}%")
        if self.data is not None and self.method_var.get() == "manual":
            self.process_image()

    def update_block_size(self, *args):
        value = self.block_size_var.get()
        # Ensure block size is odd
        if value % 2 == 0:
            value += 1
            self.block_size_var.set(value)
        self.block_size_label.config(text=f"Block Size: {value}")
        if self.data is not None and self.method_var.get() == "adaptive":
            self.process_image()

    def update_c_value(self, *args):
        value = self.c_value_var.get()
        self.c_value_label.config(text=f"C Value: {value}")
        if self.data is not None and self.method_var.get() == "adaptive":
            self.process_image()

    def update_percentile(self, *args):
        value = self.percentile_var.get()
        self.percentile_label.config(text=f"Percentile: {value:.1f}%")
        if self.data is not None and self.method_var.get() == "percentile":
            self.process_image()

    def update_log_factor(self, *args):
        value = self.log_factor_var.get()
        self.log_factor_label.config(text=f"Factor: {value:.1f}")
        if self.data is not None and self.method_var.get() == "log":
            self.process_image()

    def update_opening(self, *args):
        if self.data is not None:
            self.process_image()

    def update_closing(self, *args):
        if self.data is not None:
            self.process_image()

    def update_max_region(self, *args):
        if self.data is not None:
            self.process_image()

    def update_max_intensity(self, *args):
        value = self.max_intensity_var.get()
        self.max_intensity_label.config(text=f"Intensity Threshold: {value:.1f}%")
        if self.data is not None and self.max_region_only.get():
            self.process_image()

    def process_image(self, roi_changed=False):
        if self.data is None:
            return

        method = self.method_var.get()

        # Cache the original data to avoid making copies for each operation
        if hasattr(self, '_data_cache') and self._data_cache is not None:
            processed_data = self._data_cache
        else:
            processed_data = self.data.copy()
            self._data_cache = processed_data

        # Apply ROI optimization - only process within ROI region
        if self.roi_coords is not None or len(self.roi_polygon_points) > 0:
            # Create a masked version for processing
            roi_only = processed_data.copy()
            roi_only[~self.roi_mask] = 0
            processed_data = roi_only

        # Apply thresholding based on selected method
        if method == "manual":
            threshold_percent = self.threshold_var.get() / 100.0
            threshold_value = np.min(processed_data[processed_data > 0] if (self.roi_coords is not None or len(self.roi_polygon_points) > 0) else processed_data) + threshold_percent * (np.max(processed_data) - np.min(processed_data[processed_data > 0] if (self.roi_coords is not None or len(self.roi_polygon_points) > 0) else processed_data))
            self.binary_mask = processed_data > threshold_value

        elif method == "otsu":
            # Create a normalized version for Otsu (using non-zero data if ROI applied)
            norm_proc_data = np.zeros_like(processed_data, dtype=np.uint8)
            if self.roi_coords is not None or len(self.roi_polygon_points) > 0:
                mask = processed_data > 0
                if np.any(mask):
                    min_val = np.min(processed_data[mask])
                    max_val = np.max(processed_data[mask])
                    if max_val > min_val:
                        norm_proc_data[mask] = ((processed_data[mask] - min_val) / (max_val - min_val) * 255).astype(np.uint8)
            else:
                norm_proc_data = ((processed_data - np.min(processed_data)) / (np.max(processed_data) - np.min(processed_data)) * 255).astype(np.uint8)

            _, binary = cv2.threshold(norm_proc_data, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
            self.binary_mask = binary > 0

        elif method == "adaptive":
            # Get parameters and ensure block size is odd
            block_size = self.block_size_var.get()
            if block_size % 2 == 0:
                block_size += 1
            c_value = self.c_value_var.get()

            # Optimize normalization for ROI
            norm_proc_data = np.zeros_like(processed_data, dtype=np.uint8)
            if self.roi_coords is not None or len(self.roi_polygon_points) > 0:
                mask = processed_data > 0
                if np.any(mask):
                    min_val = np.min(processed_data[mask])
                    max_val = np.max(processed_data[mask])
                    if max_val > min_val:
                        norm_proc_data[mask] = ((processed_data[mask] - min_val) / (max_val - min_val) * 255).astype(np.uint8)
            else:
                norm_proc_data = ((processed_data - np.min(processed_data)) / (np.max(processed_data) - np.min(processed_data)) * 255).astype(np.uint8)

            binary = cv2.adaptiveThreshold(
                norm_proc_data, 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
                cv2.THRESH_BINARY, block_size, c_value
            )
            self.binary_mask = binary > 0

        elif method == "percentile":
            # Optimize percentile calculation for ROI
            percentile = self.percentile_var.get()
            if self.roi_coords is not None or len(self.roi_polygon_points) > 0:
                non_zero = processed_data[processed_data > 0]
                if non_zero.size > 0:
                    threshold = np.percentile(non_zero, percentile)
                else:
                    threshold = 0
            else:
                threshold = np.percentile(processed_data, percentile)
            self.binary_mask = processed_data > threshold

        elif method == "log":
            # Optimize log calculation for ROI
            log_data = np.log1p(processed_data)
            if self.roi_coords is not None or len(self.roi_polygon_points) > 0:
                non_zero = log_data[log_data > 0]
                if non_zero.size > 0:
                    mean_val = np.mean(non_zero)
                    std_val = np.std(non_zero)
                else:
                    mean_val = 0
                    std_val = 1
            else:
                mean_val = np.mean(log_data)
                std_val = np.std(log_data)

            factor = self.log_factor_var.get()
            threshold = mean_val + factor * std_val
            self.binary_mask = log_data > threshold

        # Apply morphological operations more efficiently
        opening_size = self.opening_var.get()
        closing_size = self.closing_var.get()

        # Only apply morphology if needed
        if opening_size > 0 or closing_size > 0:
            # If ROI is active, only apply morphology to ROI region to improve performance
            if self.roi_coords is not None or len(self.roi_polygon_points) > 0:
                # Find bounding box of ROI
                roi_indices = np.where(self.roi_mask)
                if len(roi_indices[0]) > 0:
                    y1, y2 = roi_indices[0].min(), roi_indices[0].max()
                    x1, x2 = roi_indices[1].min(), roi_indices[1].max()
                    
                    roi_binary = self.binary_mask[y1:y2+1, x1:x2+1].copy()

                    if opening_size > 0:
                        roi_binary = morphology.binary_opening(roi_binary, morphology.disk(opening_size))
                    if closing_size > 0:
                        roi_binary = morphology.binary_closing(roi_binary, morphology.disk(closing_size))

                    # Update only the ROI part of the full mask
                    self.binary_mask[y1:y2+1, x1:x2+1] = roi_binary
            else:
                # Apply to full image
                if opening_size > 0:
                    self.binary_mask = morphology.binary_opening(
                        self.binary_mask, morphology.disk(opening_size)
                    )
                if closing_size > 0:
                    self.binary_mask = morphology.binary_closing(
                        self.binary_mask, morphology.disk(closing_size)
                    )

        # Apply ROI mask to the binary mask
        if self.roi_coords is not None or len(self.roi_polygon_points) > 0:
            self.binary_mask = np.logical_and(self.binary_mask, self.roi_mask)

        # Optimize merge regions - only process if we have active regions
        if self.merge_regions.get() and np.any(self.binary_mask):
            self.merge_nearby_regions()

        # Max region filtering - only process if enabled and we have regions
        if self.max_region_only.get() and np.any(self.binary_mask):
            self.filter_max_intensity_region()

        # Optimize edge detection for better performance
        if self.roi_coords is not None or len(self.roi_polygon_points) > 0:
            # Find bounding box of ROI
            roi_indices = np.where(self.roi_mask)
            if len(roi_indices[0]) > 0:
                y1, y2 = roi_indices[0].min(), roi_indices[0].max()
                x1, x2 = roi_indices[1].min(), roi_indices[1].max()
                
                # Only compute edges in the ROI region
                roi_data = self.data[y1:y2+1, x1:x2+1].copy()
                if np.any(roi_data > 0):
                    norm_roi = (roi_data - np.min(roi_data)) / (np.max(roi_data) - np.min(roi_data))
                    roi_edges = feature.canny(norm_roi, sigma=2)

                    # Place ROI edges in full image
                    self.edges = np.zeros_like(self.data, dtype=bool)
                    self.edges[y1:y2+1, x1:x2+1] = roi_edges
                    # Apply ROI mask
                    self.edges = np.logical_and(self.edges, self.roi_mask)
                else:
                    self.edges = np.zeros_like(self.data, dtype=bool)
            else:
                self.edges = np.zeros_like(self.data, dtype=bool)
        else:
            # Process full image
            norm_data = ((self.data - np.min(self.data)) / (np.max(self.data) - np.min(self.data)))
            self.edges = feature.canny(norm_data, sigma=2)

        # Create Extracted jet efficiently
        self.Extracted_jet = np.zeros_like(self.data, dtype=float)
        self.Extracted_jet[self.binary_mask] = self.data[self.binary_mask]

        # Update visualization
        self.update_visualization()
        
    def merge_nearby_regions(self):
        if not np.any(self.binary_mask):
            return

        # Get connected components
        labeled_mask, num_features = ndimage.label(self.binary_mask)
        if num_features <= 1:  # No need to merge if there's only one region
            return

        # Get minimum size threshold
        min_size = self.merge_size_threshold.get()

        # For force merge all, use the simplest approach possible
        if self.force_merge_all.get():
            # Find all regions
            regions = measure.regionprops(labeled_mask)
            
            # Create a mask that contains only regions above the minimum size
            merged_mask = np.zeros_like(self.binary_mask, dtype=bool)
            
            # Count valid regions for potential convex hull approach
            valid_regions = []
            
            for region in regions:
                if region.area >= min_size:
                    region_mask = labeled_mask == region.label
                    merged_mask = np.logical_or(merged_mask, region_mask)
                    valid_regions.append(region)
            
            # Use convex hull if the connection strength is high enough
            if self.merge_strength.get() > 0.1 and len(valid_regions) > 1:
                # Get all points from the merged mask
                y, x = np.where(merged_mask)
                points = np.column_stack((x, y))
                
                # Only compute convex hull if we have enough points
                if len(points) > 3:
                    try:
                        from scipy.spatial import ConvexHull
                        hull = ConvexHull(points)
                        
                        # Create a polygon from the hull vertices
                        hull_vertices = points[hull.vertices]
                        
                        # Create an image filled with the hull polygon
                        hull_mask = np.zeros_like(self.binary_mask, dtype=bool)
                        rr, cc = draw_polygon(hull_vertices[:, 1], hull_vertices[:, 0], self.binary_mask.shape)
                        hull_mask[rr, cc] = True
                        
                        # Combine with original mask to ensure we don't lose valid regions
                        self.binary_mask = np.logical_or(hull_mask, merged_mask)
                        return
                    except Exception as e:
                        # Fall back to just using the merged mask if convex hull fails
                        pass
            
            # If convex hull wasn't used, just return the merged mask
            self.binary_mask = merged_mask
        
    def filter_max_intensity_region(self):
        """Keep only the region with maximum intensity - optimized version with ROI support"""
        if not np.any(self.binary_mask):
            return

        # Label connected regions
        labeled_mask, num_features = ndimage.label(self.binary_mask)
        if num_features == 0:
            return

        # Get max and threshold from ROI region only
        # If ROI is active, calculate max only within ROI, otherwise use full image
        if self.roi_coords is not None or len(self.roi_polygon_points) > 0:
            # Only consider data within ROI mask
            roi_data = self.data[self.roi_mask]
            if roi_data.size > 0:
                global_max = np.max(roi_data)
            else:
                global_max = np.max(self.data)
        else:
            global_max = np.max(self.data)
        
        threshold = self.max_intensity_var.get() / 100.0 * global_max

        # Find regions with max intensity above threshold
        valid_mask = np.zeros_like(self.binary_mask, dtype=bool)

        # More efficient iteration over regions
        for i in range(1, num_features + 1):
            region_mask = labeled_mask == i
            region_max = np.max(self.data[region_mask])
            if region_max >= threshold:
                valid_mask[region_mask] = True

        self.binary_mask = valid_mask

    def update_visualization(self):
        # Initialize LogNorm if not already created
        if self.norm_log is None:
            self.norm_log = LogNorm()

        # Clear axes without recreating the figure
        for ax in self.axes.flat:
            ax.clear()
            ax.set_xticks([])
            ax.set_yticks([])
        cmap=plt.get_cmap('solar orbiterhri_euv174')
        
        # Plot data efficiently with cached norm objects
        # Original data
        im1 = self.axes[0, 0].imshow(self.data, cmap=cmap, norm=self.norm_log, origin='lower')
        self.axes[0, 0].set_title("Original Data")

        # Manage colorbars efficiently
        if hasattr(self, 'cbar1') and self.cbar1 is not None:
            self.cbar1.update_normal(im1)
        else:
            self.cbar1 = self.fig.colorbar(im1, ax=self.axes[0, 0], fraction=0.046, pad=0.04)

        # Draw ROI rectangle or polygon
        if self.roi_coords is not None:
            x1, y1, x2, y2 = self.roi_coords
            rect = plt.Rectangle((x1, y1), x2-x1, y2-y1,
                                linewidth=2, edgecolor='lime', facecolor='none')
            self.axes[0, 0].add_patch(rect)
        elif len(self.roi_polygon_points) > 0:
            # Draw polygon
            vertices = np.array(self.roi_polygon_points)
            from matplotlib.patches import Polygon
            poly = Polygon(vertices, linewidth=2, edgecolor='lime', facecolor='none', closed=True)
            self.axes[0, 0].add_patch(poly)

        # Binary mask (efficient)
        self.axes[0, 1].imshow(self.binary_mask, cmap='gray', origin='lower')
        self.axes[0, 1].set_title("Binary Mask")

        # Extracted jet
        masked_extracted = np.ma.masked_where(~self.binary_mask, self.data)
        im3 = self.axes[1, 0].imshow(masked_extracted, cmap=cmap, norm=self.norm_log, origin='lower')
        self.axes[1, 0].set_title("Extracted Jet")
        if hasattr(self, 'cbar3') and self.cbar3 is not None:
            self.cbar3.update_normal(im3)
        else:
            self.cbar3 = self.fig.colorbar(im3, ax=self.axes[1, 0], fraction=0.046, pad=0.04)

        # Edge detection
        self.axes[1, 1].imshow(self.edges, cmap='gray', origin='lower')
        self.axes[1, 1].set_title("Edge Detection")

        # Update canvas - use draw_idle for better performance
        self.fig.tight_layout(pad=1.0)
        self.canvas.draw_idle()
        
    def save_results(self):
        if self.data is None or self.binary_mask is None:
            messagebox.showwarning("Warning", "No data to save")
            return

        # Select directory to save results
        save_dir = filedialog.askdirectory(title="Select Directory to Save Results")

        if save_dir:
            try:
                base_name = os.path.splitext(os.path.basename(self.fits_file))[0]
                output_prefix = os.path.join(save_dir, f"{base_name}_processed")

                # Create and save SunPy maps

                # Save binary mask as FITS
                mask_map = sunpy.map.Map(self.binary_mask.astype(np.uint8), self.header)
                mask_map.save(f"{output_prefix}_mask.fits", overwrite=True)

                # Save Extracted jet as FITS
                Extracted_map = sunpy.map.Map(self.Extracted_jet, self.header)
                Extracted_map.save(f"{output_prefix}_Extracted.fits", overwrite=True)

                # Save edge detection as FITS
                edges_map = sunpy.map.Map(self.edges.astype(np.uint8), self.header)
                edges_map.save(f"{output_prefix}_edges.fits", overwrite=True)

                # Save ROI information if applicable
                if self.roi_coords is not None or len(self.roi_polygon_points) > 0:
                    roi_mask = np.zeros_like(self.data, dtype=np.uint8)
                    roi_mask[self.roi_mask] = 1
                    roi_map = sunpy.map.Map(roi_mask, self.header)
                    roi_map.save(f"{output_prefix}_roi_mask.fits", overwrite=True)

                    # Also save ROI coordinates to a text file
                    with open(f"{output_prefix}_roi_coords.txt", 'w') as f:
                        if self.roi_coords is not None:
                            f.write(f"ROI Type: Rectangle\n")
                            f.write(f"ROI Coordinates: {self.roi_coords}\n")
                        else:
                            f.write(f"ROI Type: Polygon\n")
                            f.write(f"Number of vertices: {len(self.roi_polygon_points)}\n")
                            for i, vertex in enumerate(self.roi_polygon_points):
                                f.write(f"Vertex {i+1}: ({vertex[0]:.2f}, {vertex[1]:.2f})\n")

                # Save region merging information if applicable
                if self.merge_regions.get():
                    with open(f"{output_prefix}_merge_params.txt", 'w') as f:
                        f.write(f"Region Merging: Enabled\n")
                        if self.force_merge_all.get():
                            f.write(f"Force Merge All: Enabled\n")
                            f.write(f"Connection Strength: {self.merge_strength.get():.3f}\n")
                            if self.merge_strength.get() > 0.1:
                                f.write(f"Convex Hull: Enabled\n")
                            else:
                                f.write(f"Convex Hull: Disabled\n")
                        else:
                            f.write(f"Force Merge All: Disabled\n")
                            f.write(f"Distance Threshold: {self.merge_distance.get()} pixels\n")
                        f.write(f"Minimum Region Size: {self.merge_size_threshold.get()} pixels\n")

                # Save visualizations as PNG
                plt.figure(figsize=(10, 8))
                plt.imshow(self.Extracted_jet, cmap='inferno', origin='lower')
                plt.colorbar(label='Extracted Intensity')
                plt.title('Extracted Jet Structure')
                
                # Draw ROI on saved figure
                if self.roi_coords is not None:
                    x1, y1, x2, y2 = self.roi_coords
                    rect = plt.Rectangle((x1, y1), x2-x1, y2-y1,
                                        linewidth=2, edgecolor='lime', facecolor='none')
                    plt.gca().add_patch(rect)
                elif len(self.roi_polygon_points) > 0:
                    vertices = np.array(self.roi_polygon_points)
                    from matplotlib.patches import Polygon
                    poly = Polygon(vertices, linewidth=2, edgecolor='lime', facecolor='none', closed=True)
                    plt.gca().add_patch(poly)
                    
                plt.savefig(f"{output_prefix}_Extracted.png", dpi=300, bbox_inches='tight')
                plt.close()

                plt.figure(figsize=(10, 8))
                plt.imshow(self.binary_mask, cmap='gray', origin='lower')
                plt.title('Binary Mask')
                plt.savefig(f"{output_prefix}_mask.png", dpi=300, bbox_inches='tight')
                plt.close()

                messagebox.showinfo("Success", f"Results saved to {save_dir}")

            except Exception as e:
                messagebox.showerror("Error", f"Failed to save results: {str(e)}")

# Run the application
if __name__ == "__main__":
    import os
    os.environ['TK_SILENCE_DEPRECATION'] = '1'
    root = tk.Tk()
    app = InteractiveThresholdingApp(root)
    root.mainloop()