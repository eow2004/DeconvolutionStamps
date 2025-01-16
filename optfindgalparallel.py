import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from os import listdir
import multiprocessing as mp

import warnings
warnings.filterwarnings("ignore")

# Define the path to the FITS files
path = r"/home/eliasw/OneDrive/Research/Deconvolution/fits"

# Load the galaxy data and FITS file ranges
obj = pd.read_csv('withinrangeCOSMOS2020.csv', sep=',', header=0)
ranges = pd.read_csv('fitsrange.txt', sep='\t', header=0)

# Preload FITS headers and data
def preload_fits_files(path):
    fits_data = {}
    for file in listdir(path):
        if file.endswith('.fits.gz'):  # Ensure it's a FITS file
            header = WCS(fits.getheader(path + '/' + file, ext=1), fix=False)  # Get WCS header for pixel to RA/Dec conversion
            mask = fits.getdata((path + '/mask/' + file).replace('.fits.gz','_mask.fits'))  # Load only the necessary data extension
            fits_data[file] = {'wcs': header, 'mask': mask}
    return fits_data

# Function to search for galaxies in preloaded FITS files for a subset of galaxies
def search_galaxies(start_index, end_index, fits_data):
    foundgals = 0
    found_objects = []
        
    end_index = min(end_index, len(obj)) #indexing fix
    #print(f"Processing indices from {start_index} to {end_index}")

    for i in range(start_index, end_index):
       # print(f"Checking index {i}, Galaxy ID {obj['Classic'][i]}") 
        for file, file_data in fits_data.items():  # Loop through preloaded FITS data
            ramin = float(ranges[ranges['File'] == file]['RAmin'])
            ramax = float(ranges[ranges['File'] == file]['RAmax'])
            decmin = float(ranges[ranges['File'] == file]['DECmin'])
            decmax = float(ranges[ranges['File'] == file]['DECmax'])
            # Check if the galaxy is within the RA/Dec range of the file
            if ramin <= obj['RAJ2000'][i]:
                if obj['RAJ2000'][i] <= ramax:
                    if decmin <= obj['DEJ2000'][i]:
                        if obj['DEJ2000'][i] <= decmax:
                            w = file_data['wcs']  # Get the preloaded WCS header
                            xpix, ypix = w.all_world2pix(obj['RAJ2000'][i], obj['DEJ2000'][i], 1)
                            # Check if the pixel coordinates are within bounds
                            if (0 <= xpix <= 9600 and 0 <= ypix <= 12455):
                                if (file_data['mask'][int(ypix),int(xpix)] == True):  # Galaxy found
                                    foundgals += 1
                                    #print(f"Galaxy found at index {i}, ID {obj['Classic'][i]} in file {file}, at coords {obj['RAJ2000'][i]} RA and {obj['DEJ2000'][i]} DEC")
                                    found_objects.append([obj['Classic'][i], obj['RAJ2000'][i], obj['DEJ2000'][i], file])  # Save RA, Dec, and file
                                    break  # Stop searching this file if galaxy found
    
    return foundgals, found_objects

# Multi-core processing: split the galaxy list into chunks
if __name__ == '__main__':
    num_cores = mp.cpu_count()  # Use all available cores
    chunk_size = len(obj) // num_cores  # Determine the size of each chunk

    # Preload all FITS files and headers before starting the search
    fits_data = preload_fits_files(path)

    # Create a pool of workers and distribute the work
    with mp.Pool(num_cores) as pool:
        results = [pool.apply_async(search_galaxies, args=(i, i + chunk_size, fits_data)) for i in range(0, len(obj), chunk_size)]
        
        found_gals = 0
        found_objects = []

        for result in results:
            chunk_found_gals, chunk_found_details = result.get()  # Get both elements from the result tuple
            found_gals += chunk_found_gals  # Sum the galaxy count
            found_objects.extend(chunk_found_details)  # Add the found details to the main list

        # Write the results to a file
        with open('found_galaxies.txt', 'w') as out:
            out.write(f'Found {found_gals} galaxies out of {len(obj)}\n')
            out.write('ID\tRAJ2000\tDEJ2000\tFile\n')  # Header for the galaxy details
            for gal in found_objects:
                out.write(f'{gal[0]}\t{gal[1]}\t{gal[2]}\t{gal[3]}\n')  # Write each galaxy's details (RA, Dec, and file)
            out.flush()
