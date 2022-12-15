# CardiacNoiseR2s

# Authors
Quentin Raynaud, 2022

Laboratory for Neuroimaging Research,
Lausanne University Hospital & University of Lausanne, Lausanne, Switzerland

Copyright (C) 2022 Laboratory for Neuroimaging Research

# Content

The scripts used in CardiacNoiseR2s are the ones used to perform the characterization of cardiac-induced noise, in the article "Characterization of cardiac-induced noise in R2* maps of the brain", from Quentin Raynaud, Giulia Di Domenicantonio, Jérome Yerly, Thomas Dardano, Ruud B. van Heeswijk and Antoine Lutti.

A subset of data used in the original publication for computations of the results can be downloaded on Zenodo DOI:10.5281/zenodo.7428605.
These data contains 5D k-space data (kx/ky/kz/cardiac phase/TE) from one of the participants. 

"CardiacNoiseR2s.zip" contain the data that need to be used with the code.

"CardiacNoiseR2s_validation.zip" contains the main results outputted by the code, and can be used for validation.

# Main functions

qr_NoiseCharacLowRes_Main is the main function that needs to be run by inputting the folder where the subject data are located.
It requires the 5D k-space data (Truek.mat), the labels (Labels.mat) as well as header information for SPM (SPMheader.mat).

This function will run the main functions in the following order: qr_CardiacNoiseR2s_Analysis, qr_CardiacNoiseR2s_R2sCharact, qr_CardiacNoiseR2s_SensitiveFreqs, and qr_CardiacNoiseR2s_ShowResults.

## qr_CardiacNoiseR2s_Analysis
First main function used to model the cardiac-induced noise in image and k-space with fundamental and first harmonic of the Fourier series.
It saves the Fourier series weights (betas) as nifties for visualization using SMP12, and as .mat to be loaded with matlab. It also saves the modeled cardiac-induced noise in k-space.

## qr_CardiacNoiseR2s_R2sCharact
Second main function used to compute maps of R2*/RMSE, R2*/RMSE standard deviation across the cardiac phase of the raw data, the raw data when the modeled cardiac noise is removed, and the modeled cardiac noise only.

## qr_CardiacNoiseR2s_SensitiveFreqs
Third main function used to find the cardiac noise sensitive k-space region. It computes maps of R2*/RMSE, R2*/RMSE standard deviation across the cardiac phase of the raw data and the modeled cardiac noise only, when removing the fundamental and first harmonic of the Fourier series on a specific region of k-space.

## qr_CardiacNoiseR2s_ShowResults
Final main function that loads all the relevant results and makes images from them for visualization. The different images are saved in \Results.

# REQUIREMENTS

The scripts were run under MATLAB 2021a.
The scripts require SPM12 to function properly (avaliable here: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/).

