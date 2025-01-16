# DeconvolutionStamps
In preparating for comparing JWST data with deconvolved images of galaxies, I have done work to prepare tools to image galaxies based on coordinates.

## Objective
### JWST Data Release (COSMOS-Web) Stamps
I am going to do research testing and evaluating a new deconvolution tool on a specific type of galaxies. To give an accurate evaluation of the images produced by the deconvolution model, we need a baseline to compare them against.

Ground based radio observations are fairly availible and Hubble images stand to be a good baseline to compare against; however, with a new state-of-the-art telescope in the form of the James-Webb Space Telescope, coupled with the fact that the deconvolution model was trained and modeled around JWST images, I figured I should play with the new public data release.

Currently, I have worked through 2 Jupyter Notebooks which prepare singular images of galaxies in multiple bands, able to be used for comparison.


## [stamps.ipynb](stamps.ipynb)
Stamps handles reading in singular .fits images of JWST MIRI and NIRCAM data and preparing "stamps" of singular galaxies in a .png, .fits, and mosaic format.

When working with the whole data release, Stamps, shows how I handled finding galaxies across all of the mosaic .fits files from COSMOS-Web. The end result is [optfindgalparallel.py](optfindgalparallel.py), a script that takes a list of target galaxies and their coordinates. The script determines whether the coordinates lies within any of the mosaics, and whether a mask image pre-prepared has data at that point. This exports a list of "found" galaxies and what mosaic it is located in, which can be used (later in Imaging) to iterate over and have a stamp made.

Stamps also holds code I made to find the borders of the mosaics (used in optfindgalparallel).

Stamps shows my budding proficiency in astro data handling and python libraries, and showcases linearly my improvement and optimization of my code. I have tried to separate the useful code from the errored code, but it may take some cleaning.

I plan to do some improvements necessary once I prepare to focus on a specific subset of galaxies:
* Add a culling parameter to optfindgalparallel to exclude coordinates with values too high (stars)
* Add a culling parameter to optfindgalparallel to exclude coordinates with values too low (galaxies that may or may not be there, but are too faint to do meaningfull analysis with)
* Add a culling parameter to optfindgalparallel to exclude galaxies too close to another.
* Etc.
I may also reformat it into a form where it can be called as a function in a code that images the galaxies, but only if it is more effecient than calling a pre-made catalog like it currenly does.


## [imaging.ipynb](imaging.ipynb)
Imaging handles work post-identifying galaxies that can be imaged. Imaging creates lists of galaxies in each mosaic in the COSMOS-Web data, and then produces stamps in the 115, 150, 277, and 444 micron bands for each galaxy. It displays these in a series next to each other. (This can and might be modified into producing the stamps as a singular RGB .png image)
![A2_649770_allbands](https://github.com/user-attachments/assets/1046ab5d-03cc-4959-b161-0e33876c6119)
![ds9](https://github.com/user-attachments/assets/29a40a29-9a37-4d50-bd32-83bb83050252)


Imaging's second act focuses on the catalog COSMOS2020 (Weaver et. al 2020), and connecting the data from such. 

All of the "target" galaxies in this project come from COSMOS2020, as the data attached to them will be instrumental in identifying the demographics.

After producing the stamps for each band, Imaging can save them in a compressed .fits directory with a fits table with data for the galaxy. This data is from the COSMOS2020 catalog. This may be edited and/or scrapped, depending on whether I use the data to make checks and cuts in Stamps instead of Imaging.

Imaging shows my work with using data handling and imaging to produce composite images of a galaxy, and the latter part showcases my abilities in handling .fits files, headers, and comments to effectively store data for astronomical use.

## Further Work
Some improvements needed before futher work, or alongside futherwork are:
* !!! Handle images are cropped by the edges of the mosaics, which results in an error and an unusable image !!!
* Adding redshift labels on the produced stamps.
* Changing axes on produced stamps to RA/DEC coordinates.
* Producing code to recursively create RGB images from the multi-band stamps.


My next objectives are to:
* Find a demographic to compare deconvolved images against.
* Cull and produce JWST stamps of these galaxies using the COSMOS2020 catalog and COSMOS-Web images.
* Use the deconvolution model to produce enhanced images of these galaxies.
* Compare and analyize the work of the model.
* Make observations of the demographic using deconvolved and JWST images.
