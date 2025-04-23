# CREDIT Aavik Wadivkar

import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog, messagebox
from astropy.io import fits
from scipy.ndimage import gaussian_filter
from skimage.transform import resize
from skimage.io import imread
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# Default Constants
DEFAULTS = {
    "Pixel Size (um)": 3.76e-6,
    "Sensor Width (px)": 9568,
    "Sensor Height (px)": 6380,
    "Aperture": 0.111,
    "QE": 0.6,
    "Wavelength (nm)": 640,
    "Dark Current (e-)": 0.56,
    "Saturation Capacity (e-)": 51000,
    "Readout Noise (e-)": 1,
    "Field of View (deg)": 10,
    "Max Magnitude": 20,
    "Min Magnitude": 12,
    "Zero Point": 18,
    "PSF (sigma)": 3,
    "Exposure Time": 10,
    "Num of Stars": 1000,
    "Trail Length (px)": 10,
    "Drift Angle (deg)": 0,
    "Cosmic Ray Count": 5,
    "Cosmic Ray Max Length": 20,
    "Cosmic Ray Intensity (e-)": 5000,
    "Sky Background Rate (e-/px/s)": 0.1,
    "Image Scale Factor": 1.0,
    "Image Magnitude": 1.0
}

## Helper Functions

def apply_binning(image, bin_size=3):
    """Bin the image in non-overlapping blocks of bin_size x bin_size pixels."""
    h, w = image.shape
    # Trim the image so that dimensions are divisible by bin_size:
    h_trim = h - (h % bin_size)
    w_trim = w - (w % bin_size)
    trimmed = image[:h_trim, :w_trim]
    # Reshape and sum over the binning blocks.
    binned = trimmed.reshape(h_trim // bin_size, bin_size, w_trim // bin_size, bin_size).sum(axis=(1,3))
    return binned

def add_cosmic_rays(image, num_rays=5, max_length=20, intensity=5000):
    ## Add simulated cosmic ray events as short bright streaks."""
    for _ in range(num_rays):
        # Choose a random starting pixel.
        start_x = np.random.randint(0, image.shape[1])
        start_y = np.random.randint(0, image.shape[0])
        # Random length and direction.
        length = np.random.randint(1, max_length)
        angle = np.random.uniform(0, 2 * np.pi)
        for i in range(length):
            x = int(start_x + i * np.cos(angle))
            y = int(start_y + i * np.sin(angle))
            if 0 <= x < image.shape[1] and 0 <= y < image.shape[0]:
                image[y, x] += intensity
    return image

def add_sky_background(image, background_rate, exposure_time):
    """Add a uniform sky background (in electrons) across the image."""
    # In a more detailed model, you’d need the sun’s elevation, atmospheric conditions, etc.
    # Check skybackcalc.py for a (presumably) more accurate calculation.
    return image + background_rate * exposure_time

imported_image = []
image_preview_canvas = None

def import_image():
    global imported_image
    filename = filedialog.askopenfilename(
        title="Select a PNG, JPEG, or FITS File",
        filetypes=[
            ("PNG Files", "*.png"),
            ("FITS Files", "*.fits"),
            ("JPG Files", "*.jpg"),
            ("JPEG Files", "*.jpeg"),
            ("All Files", "*.*"),
        ]
    )

    if not filename:  # If user cancels, exit function safely
        return

    try:
        if filename.lower().endswith(".fits"):
            with fits.open(filename) as hdul:
                imported_image = hdul[0].data.astype(float)
        else:
            imported_image = imread(filename, as_gray=True).astype(float)

        if imported_image is None or imported_image.size == 0:
            raise ValueError("Loaded image is empty or invalid.")

        messagebox.showinfo("Success", "Image imported successfully!")
        update_image_preview()

    except Exception as e:
        messagebox.showerror("Error", f"Failed to load image: {str(e)}")
        imported_image = None  # Reset in case of failure

def update_image_preview():
    global image_preview_canvas, imported_image
    if imported_image is not None:
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.imshow(imported_image, cmap='gray', origin='lower')
        ax.set_title("Imported Image Preview")
        ax.axis('off')
        
        if image_preview_canvas:
            image_preview_canvas.get_tk_widget().destroy()
        
        image_preview_canvas = FigureCanvasTkAgg(fig, master=preview_frame)
        image_preview_canvas.get_tk_widget().pack()
        image_preview_canvas.draw()

def calculate_image_flux(image_magnitude, sensor_params):
    """Calculate appropriate photon flux per pixel for the imported image."""
    total_flux = 10 ** (-0.4 * (image_magnitude - sensor_params["Zero Point"]))
    total_photons = total_flux * sensor_params["Exposure Time"]  # Adjust based on exposure time
    total_photons *= sensor_params["QE"]  # Adjust based on quantum efficiency
    # total_photons *= 0.1 # Adjust based on \delta 30 nm bandpass filter
    pixel_flux = total_photons / np.sum(imported_image)  # Normalize across image pixels
    return pixel_flux

def add_imported_image(image, scale_factor=1.0, image_magnitude=10.0, sensor_params=None):
    global imported_image
    if imported_image is not None:
        original_h, original_w = imported_image.shape
        sensor_h, sensor_w = image.shape

        # Compute new dimensions while maintaining aspect ratio
        new_width = int(sensor_w * scale_factor)
        new_height = int(original_h * (new_width / original_w))

        # Ensure the resized image does not exceed sensor dimensions
        new_height = min(new_height, sensor_h)
        new_width = min(new_width, sensor_w)

        resized_image = resize(imported_image, (new_height, new_width), anti_aliasing=True)
        pixel_flux = calculate_image_flux(image_magnitude, sensor_params)
        resized_image *= pixel_flux  # Scale the image based on calculated flux

        # Compute centering positions
        start_x = max(0, (sensor_w - new_width) // 2)
        start_y = max(0, (sensor_h - new_height) // 2)

        # Ensure the final cropped region does not exceed bounds
        end_x = start_x + new_width
        end_y = start_y + new_height

        # Crop to ensure compatibility with NumPy broadcasting
        image[start_y:end_y, start_x:end_x] += resized_image[:end_y - start_y, :end_x - start_x]

    return image

### PSF addition ###

def add_psf(image, x, y, flux, sigma):
    """Add a 2D Gaussian PSF to the image at (x, y) with the given flux."""
    size = int(6 * sigma)
    # Create grid centered on 0
    y_indices, x_indices = np.meshgrid(np.arange(-size//2, size//2+1),
                                       np.arange(-size//2, size//2+1), indexing='ij')
    psf = np.exp(-(x_indices**2 + y_indices**2) / (2 * sigma**2))
    psf /= psf.sum()
    ix, iy = int(y), int(x)
    # Determine sub-image boundaries
    if 0 <= ix < image.shape[0] and 0 <= iy < image.shape[1]:
        x_start, x_end = max(0, ix - size//2), min(image.shape[0], ix + size//2+1)
        y_start, y_end = max(0, iy - size//2), min(image.shape[1], iy + size//2+1)
        sub_psf = psf[:x_end-x_start, :y_end-y_start]
        image[x_start:x_end, y_start:y_end] += flux * sub_psf

### Main image generation function ###

def generate_image(params, binning=False, cosmic_rays=False, sky_background=False, moving_exposures=False, snr_calc=False, import_image=False):
    sensor_h = int(params["Sensor Height (px)"])
    sensor_w = int(params["Sensor Width (px)"])
    image = np.zeros((sensor_h, sensor_w))
    
    # Generate star positions and magnitudes.
    num_stars = int(params["Num of Stars"])
    x_positions = np.random.uniform(0, sensor_w, num_stars)
    y_positions = np.random.uniform(0, sensor_h, num_stars)
    magnitudes = np.random.uniform(params["Min Magnitude"], params["Max Magnitude"], num_stars)
    
    if moving_exposures:
        # Divide the exposure into subexposures to simulate camera drift.
        num_steps = 10  # Can be a parameter?

        trail_length = params["Trail Length (px)"]
        drift_angle_rad = np.deg2rad(params["Drift Angle (deg)"])
        dx = trail_length * np.cos(drift_angle_rad) / (num_steps - 1)
        dy = trail_length * np.sin(drift_angle_rad) / (num_steps - 1)
        for step in range(num_steps):
            sub_exposure_time = params["Exposure Time"] / num_steps
            # For each subexposure, add the star signals with a small positional offset
            for x, y, mag in zip(x_positions, y_positions, magnitudes):
                flux = 10 ** (-0.4 * (mag - params["Zero Point"]))
                photons = flux * sub_exposure_time
                electrons = np.random.poisson(photons * params["QE"])
                add_psf(image, x + dx * step, y + dy * step, electrons, params["PSF (sigma)"])
    else:
        # Normal (static) exposure
        for x, y, mag in zip(x_positions, y_positions, magnitudes):
            flux = 10 ** (-0.4 * (mag - params["Zero Point"]))
            photons = flux * params["Exposure Time"]
            electrons = np.random.poisson(photons * params["QE"])
            add_psf(image, x, y, electrons, params["PSF (sigma)"])
    
    # Add imported image if available
    if import_image:
        image = add_imported_image(image, params["Image Scale Factor"], params["Image Magnitude"], params)

    # Add cosmic rays if toggled
    if cosmic_rays:
        image = add_cosmic_rays(image,
                                 num_rays=int(params["Cosmic Ray Count"]),
                                 max_length=int(params["Cosmic Ray Max Length"]),
                                 intensity=int(params["Cosmic Ray Intensity (e-)"]))
    
    # Add sky background if toggled
    if sky_background:
        image = add_sky_background(image,
                                   background_rate=params["Sky Background Rate (e-/px/s)"],
                                   exposure_time=params["Exposure Time"])
    
    # Add dark noise and readout noise
    dark_noise = np.random.poisson(params["Dark Current (e-)"] * params["Exposure Time"], image.shape)
    readout_noise = np.random.normal(params["Readout Noise (e-)"], 1.5, image.shape).astype(int)
    image += dark_noise + readout_noise
    
    # Clip the image to sensor's saturation capacity
    image = np.clip(image, 0, params["Saturation Capacity (e-)"]).astype(int)
    
    # Apply 3x3 binning if toggled
    if binning:
        print('Applying 3x3 binning')
        image = apply_binning(image, bin_size=3)
    
    if snr_calc:
        print('Calculating SNR')
        # Calculate the signal-to-noise ratio.
        
        signal = image.mean()
        noise = np.std(image)
        snr = 10*np.log10(signal / noise)
        print(f"Signal: {signal}, Noise: {noise}, SNR: {snr}")

        # This is a simple estimate assuming Poisson noise for CCD/CMOS.
        signal = params["Exposure Time"] * params["QE"] * 10**(-0.4 * (params["Min Magnitude"] - params["Zero Point"]))
        noise = np.sqrt(signal + int(sky_background_var.get()) * params["Sky Background Rate (e-/px/s)"] + params["Dark Current (e-)"] * params["Exposure Time"] + params["Readout Noise (e-)"]**2)
        snr = signal / noise
        print("Expected Signal: ", signal, "Expected Noise: ", noise, "Expected SNR: ", snr)

    return image

### Image saving ###

def save_image(image, filename, format):
    if format == 'png':
        plt.imsave(filename, image, cmap='gray')
    elif format == 'fits':
        fits.writeto(filename, image, overwrite=True)

### GUI Control ###

def run_simulation():
    # Retrieve parameters from the text entries
    params = {key: float(entries[key].get()) if entries[key].get() else DEFAULTS[key]
              for key in DEFAULTS}
    image = generate_image(params,
                           binning=binning_var.get(),
                           cosmic_rays=cosmic_rays_var.get(),
                           sky_background=sky_background_var.get(),
                           moving_exposures=moving_exposures_var.get(),
                           snr_calc=snr_calc_var.get(),
                           import_image=import_image_var.get())
    plt.imshow(image, cmap='gray', origin='lower')
    plt.colorbar(label='Electron Count')
    plt.title('Simulated CMOS Image')
    plt.show()    

    plt.hist(image.flatten(), bins=min(image.max(), 100), range=(0, image.max()))
    plt.yscale('log')
    plt.show()

def save_file():
    filename = filedialog.asksaveasfilename(defaultextension=".png",
                                            filetypes=[("PNG Image", "*.png"), ("FITS File", "*.fits")])
    if filename:
        format = "fits" if filename.endswith(".fits") else "png"
        params = {key: float(entries[key].get()) if entries[key].get() else DEFAULTS[key]
                  for key in DEFAULTS}
        image = generate_image(params,
                               binning=binning_var.get(),
                               cosmic_rays=cosmic_rays_var.get(),
                               sky_background=sky_background_var.get(),
                               moving_exposures=moving_exposures_var.get(),
                               snr_calc = snr_calc_var.get(),
                               import_image = import_image_var.get())
        save_image(image, filename, format)
        messagebox.showinfo("Success", f"Image saved as {filename}")

### GUI Setup ###

root = tk.Tk()
root.title("CMOS Image Simulation GUI")

# Create parameter entries frame
param_frame = tk.Frame(root)
param_frame.grid(row=0, column=0, padx=10, pady=10)

entries = {}
# List default parameters in GUI
for i, (key, value) in enumerate(DEFAULTS.items()):
    tk.Label(param_frame, text=key).grid(row=i, column=0, sticky="e")
    entries[key] = tk.Entry(param_frame, width=12)
    entries[key].grid(row=i, column=1)
    entries[key].insert(0, str(value))

# Create frame for feature toggles
toggle_frame = tk.LabelFrame(root, text="Simulation Features", padx=10, pady=10)
toggle_frame.grid(row=0, column=1, padx=10, pady=10, sticky="n")

binning_var = tk.BooleanVar(value=False)
cosmic_rays_var = tk.BooleanVar(value=False)
sky_background_var = tk.BooleanVar(value=False)
moving_exposures_var = tk.BooleanVar(value=False)
snr_calc_var = tk.BooleanVar(value=False)
import_image_var = tk.BooleanVar(value=False)

tk.Checkbutton(toggle_frame, text="3x3 Binning", variable=binning_var).pack(anchor="w")
tk.Checkbutton(toggle_frame, text="Add Cosmic Rays", variable=cosmic_rays_var).pack(anchor="w")
tk.Checkbutton(toggle_frame, text="Add Sky Background", variable=sky_background_var).pack(anchor="w")
tk.Checkbutton(toggle_frame, text="Simulate Moving Exposures", variable=moving_exposures_var).pack(anchor="w")
tk.Checkbutton(toggle_frame, text="Calculate SNR", variable=snr_calc_var).pack(anchor="w")
tk.Checkbutton(toggle_frame, text="Import Image", variable=import_image_var).pack(anchor="w")

# Create frame for action buttons
button_frame = tk.Frame(root)
button_frame.grid(row=1, column=0, columnspan=2, pady=10)

tk.Button(button_frame, text="Run Simulation", command=run_simulation).grid(row=0, column=0, padx=5)
tk.Button(button_frame, text="Save Image", command=save_file).grid(row=0, column=2, padx=5)
tk.Button(button_frame, text="Import Image", command=import_image).grid(row=0, column=3, padx=5)
preview_frame = tk.LabelFrame(root, text="Imported Image Preview", padx=10, pady=10)
preview_frame.grid(row=0, column=2, padx=10, pady=10, sticky="n")

root.mainloop()