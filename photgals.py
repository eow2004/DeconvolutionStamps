import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from os import listdir
import math
import gc
import multiprocessing as mp
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# Define the path to the FITS files
path = r"/home/eliasw/OneDrive/Research/Deconvolution/fits/miri"

obj = pd.read_csv('input/asu.tsv', sep='\t', header=0,skiprows=51)

# data handling 
obj = obj.dropna() #drop all rows with missing
#obj = obj.dropna(subset=['FUVISTAH', 'e_FUVISTAH', 'FUVISTAY', 'e_FUVISTAY', 'FUVISTAKs', 'e_FUVISTAKs',]) # drop rows with blank columns

obj = obj.apply(pd.to_numeric, errors='coerce') # convert all rows
#cols_to_convert = ['FUVISTAH', 'e_FUVISTAH', 'FUVISTAY', 'e_FUVISTAY', 'FUVISTAKs', 'e_FUVISTAKs']
#obj[cols_to_convert] = obj[cols_to_convert].apply(pd.to_numeric, errors='coerce') #convert to numeric, specifc columns

obj['SNR_FUVISTAH']=(obj['FUVISTAH']/obj['e_FUVISTAH'])
obj['SNR_FUVISTAY']=(obj['FUVISTAY']/obj['e_FUVISTAY'])
obj['SNR_FUVISTAKs']=(obj['FUVISTAKs']/obj['e_FUVISTAKs'])

magFlag = (((25.5 > obj['UVISTAHmagAuto']) & (obj['UVISTAHmagAuto'] > 24.5)) & ((25.5 > obj['UVISTAHmagAuto']) & (obj['UVISTAHmagAuto'] > 24.5)) & ((25.5 > obj['UVISTAKsmagAuto']) & (obj['UVISTAKsmagAuto'] > 24.5)))
sizeFlag = ((obj['radFlux'] > 3) & (6.5 > obj['radFlux']))
snrFlag = ((obj['SNR_FUVISTAH'] > 3) & (obj['SNR_FUVISTAKs'] > 3) & (obj['SNR_FUVISTAY'] > 3))

obj = obj[magFlag & sizeFlag & snrFlag == True].reset_index()
ranges = pd.read_csv('fits/fitsrange.txt', sep=',', header=0)

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
            ramin = (ranges[ranges['File'] == file]['RAmin'])
            ramax = (ranges[ranges['File'] == file]['RAmax'])
            decmin = (ranges[ranges['File'] == file]['DECmin'])
            decmax = (ranges[ranges['File'] == file]['DECmax'])
            ramin = float(ramin)
            ramax = float(ramax)
            decmin = float(decmin)
            decmax = float(decmax)
            # Check if the galaxy is within the RA/Dec range of the file
            if ramin <= obj['RAJ2000'][i]:
                if obj['RAJ2000'][i] <= ramax:
                    if decmin <= obj['DEJ2000'][i]:
                        if obj['DEJ2000'][i] <= decmax:
                            w = file_data['wcs']  # Get the preloaded WCS header
                            xpix, ypix = w.all_world2pix(obj['RAJ2000'][i], obj['DEJ2000'][i], 1)
                            # Check if the pixel coordinates are within bounds
                            if (0 <= xpix <= 9600 and 0 <= ypix <= 12455):
                                if (file_data['mask'][int(ypix),int(xpix)] == True):  # if data at position in mask
                                    foundgals += 1
                                    #print(f"Galaxy found at index {i}, ID {obj['Classic'][i]} in file {file}, at coords {obj['RAJ2000'][i]} RA and {obj['DEJ2000'][i]} DEC")
                                    found_objects.append([obj['Classic'][i], obj['RAJ2000'][i], obj['DEJ2000'][i],obj['EZzphot'][i], obj['radFlux'][i], obj['SNR_FUVISTAH'][i], obj['UVISTAHmagAuto'][i], obj['e_UVISTAHmagAuto'][i], obj['FUVISTAH'][i], obj['e_FUVISTAH'][i], obj['SNR_FUVISTAKs'][i], obj['UVISTAKsmagAuto'][i], obj['e_UVISTAKsmagAuto'][i], obj['FUVISTAKs'][i], obj['e_FUVISTAKs'][i], obj['SNR_FUVISTAY'][i], obj['UVISTAHmagAuto'][i], obj['e_UVISTAHmagAuto'][i], obj['FUVISTAY'][i], obj['e_FUVISTAY'][i], file])  # Save RA, Dec, and file
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
        with open('output/crit_found_galaxies.txt', 'w') as out:
            #out.write(f'Found {found_gals} galaxies out of {len(obj)}\n')
            out.write('Classic,RAJ2000,DEJ2000,EZzphot,radFlux,SNR_FUVISTAH,UVISTAHmagAuto,e_UVISTAHmagAuto,FUVISTAH,e_FUVISTAH,SNR_FUVISTAKs,UVISTAKsmagAuto,e_UVISTAKsmagAuto,FUVISTAKs,e_FUVISTAKs,SNR_FUVISTAY,UVISTAHmagAuto,e_UVISTAHmagAuto,FUVISTAY,e_FUVISTAY,Section\n')  # Header for the galaxy details
            for gal in found_objects:
                out.write(f'{gal[0]},{gal[1]},{gal[2]},{gal[3]},{gal[4]},{gal[5]},{gal[6]},{gal[7]},{gal[8]},{gal[9]},{gal[10]},{gal[11]},{gal[12]},{gal[13]},{gal[14]},{gal[15]},{gal[16]},{gal[17]},{gal[18]},{gal[19]},{gal[20].replace('mosaic_miri_f770w_COSMOS-Web_60mas_',"").replace('_v0_5_i2d.fits.gz','')}\n')  # Write each galaxy's details (ID,RA, Dec, redshift, fluxradius, photometrics, and file)
            out.flush()

print('Done with gal search\nClearing fits memory.')

# Get all open file descriptors
all_open_files = [obj for obj in gc.get_objects() if isinstance(obj, fits.hdu.hdulist.HDUList)]

# Close all open FITS files
for open_file in all_open_files:
    open_file.close()

gals = pd.read_csv('output/crit_found_galaxies.txt')

print('Beginning gal image')

def decimal_2_dec(decimal):
    sign = "+" if decimal >= 0 else "-"
    degrees = math.floor(abs(decimal))
    minutes = math.floor((abs(decimal) - degrees) * 60)
    seconds = ((abs(decimal) - degrees) * 60 - minutes) * 60
    return f"{sign}{int(degrees)}°{int(minutes)}'{seconds:.2f}"

def decimal_2_ra(decimal):
    hours = int(decimal / 15)  # 360 degrees / 24 hours = 15 degrees/hour
    minutes = int((decimal / 15 - hours) * 60)
    seconds = (decimal / 15 - hours - minutes / 60) * 3600

    return f"{hours:02}:{minutes:02}:{seconds:06.2f}"

u=0

# read in section lists
sections = {}
for sec in ['A2', 'A6', 'E1', 'A10', 'A1', 'A9', 'A4', 'A8', 'A7', 'A3', 'A5']:
    sections.update({sec : (gals[gals['Section']==sec]).reset_index()})

for sec in ['A2', 'A6', 'E1', 'A10', 'A1', 'A9', 'A4', 'A8', 'A7', 'A3', 'A5']:
    pics = {}
    for band in ['115','150','277']:
        pics.update({band : fits.getdata('fits/nircam/'+sec+'/'+'mosaic_nircam_f'+band+'w_COSMOS-Web_60mas_'+sec+'_v0_5_i2d.fits.gz',ext=1)})
    w = WCS(fits.getheader('fits/nircam/'+sec+"/mosaic_nircam_f150w_COSMOS-Web_60mas_"+sec+'_v0_5_i2d.fits.gz',ext=1))
    for i in range(len(sections[sec])):
        num_images  = 3
        num_cols    = 3
        num_rows    = 1
        # Create a grid of subplots.
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(17,6))
        list_axes = list(axes.flat)
        for n,band in enumerate(['115','150','277']): # for each object
            xpix, ypix = w.all_world2pix(sections[sec]['RAJ2000'][i],sections[sec]['DEJ2000'][i],1) # convert from RA/DEC to pixel position
            xpix = round(float(xpix))
            ypix = round(float(ypix))
            rpix = round(float(25)) # COSMOS pix * pix/deg = deg / (pix/deg) = pix radius
            #print(xpix,ypix,rpix) # debug print file values
            #print(f'ra:{sections[sec]['RAJ2000'][i]},dec:{sections[sec]['DEJ2000'][i]}') # debug: print object ra/dec
            #print(f'ra:{decimal_2_ra(sections[sec]['RAJ2000'][i])},dec:{decimal_2_dec(sections[sec]['DEJ2000'][i])}') # debug: print object ra/dec in formatted string
            # get image bounds for stamp
            left = int((xpix-rpix))
            right = int((xpix+rpix))
            down = int((ypix-rpix))
            up = int((ypix+rpix))
            #print(left,right,up,down) # debug print image bounds
            # get coordinate range from left to right and up to down on the image
            leftra, downdec = w.all_pix2world(left,down,1)
            rightra, updec = w.all_pix2world(right,up,1)
            # print(leftra,rightra,updec,downdec) # debug: print image bounds coordinates
            # crop image with matplotlib
            list_axes[n].imshow(pics[band][down:up,left:right],cmap='gray',vmin=0,vmax=0.725,origin='lower',extent=[0,1,0,1],aspect='auto')
            list_axes[n].set_aspect('equal', adjustable='box')
            # print(f'xlabels:{np.linspace(rightra,leftra,3)}') # debug, print calculated image bound coordinates spaced
            list_axes[n].set_xticks([0,0.5,1],labels=[decimal_2_ra(np.linspace(leftra,rightra,3)[i]) for i in range(3)])
            list_axes[n].set_yticks([0,0.5,1],labels=[decimal_2_dec(np.linspace(downdec,updec,3)[i]) for i in range(3)])
            list_axes[n].set_title(f'{sec}_{sections[sec]['Classic'][i]}_f{band}w')
        list_axes[1].text(0.5, -0.20, 
            f"H band (mag): {round(sections[sec]['UVISTAHmagAuto'][i],5)}    Ks band (mag): {round(sections[sec]['UVISTAKsmagAuto'][i],5)}    Y band (mag): {round(sections[sec]['UVISTAHmagAuto'][i],5)}\nH band SNR: {round(sections[sec]['SNR_FUVISTAH'][i],5)}    Ks band SNR: {round(sections[sec]['SNR_FUVISTAKs'][i],5)}    Y band SNR: {round(sections[sec]['SNR_FUVISTAY'][i],5)}\nR50: {round(sections[sec]['radFlux'][i],2)}    z: {round(sections[sec]['EZzphot'][i],3)}", 
            horizontalalignment='center', wrap=True, bbox = (dict(facecolor='white', alpha=1)))
        plt.tight_layout()
        plt.savefig(f'output/Comparisons/{sec}_{sections[sec]['Classic'][i]}_allbands',dpi=75)
        print(u)
        u+=1